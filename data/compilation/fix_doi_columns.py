#!/usr/bin/env python3
"""
Fix the leaf wax calibration CSV:
- Drop data_source and audit_flags columns
- Replace Ladd compilation DOIs with original reference DOIs (Crossref lookup)
- Look up missing DOIs for the User compilation rows (Crossref lookup)
- Save back to the same file

Critical rule: NEVER fabricate a DOI. Only use a DOI that Crossref returns
and that plausibly matches the expected paper based on author + year + topic.
"""

from __future__ import annotations

import csv
import json
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

CSV_PATH = Path(
    "/Users/abradley/Desktop/proxy_uncertainty/_leafwax_paper/_for_GCA/global_data_c29_final.csv"
)
LADD_COMPILATION_DOI = "10.1029/2020JG005891"
LADD_PAPER_SOURCE = "Ladd et al. 2021"  # This one keeps the Ladd DOI

KEPT_COLS = [
    "source",
    "compilation",
    "location",
    "latitude",
    "longitude",
    "elevation",
    "d2H_precip",
    "d2H_precip_err",
    "d2H_precip_year",
    "chain",
    "d2H_wax",
    "d2H_wax_err",
    "DOI",
    "sample_type",
]

CROSSREF_URL = "https://api.crossref.org/works"
USER_AGENT = "leafwax-doi-fixer/1.0 (mailto:abradley@wustl.edu)"

# Topical terms that suggest a leaf wax / n-alkane / hydrogen isotope paper.
TOPICAL_KEYWORDS = [
    "leaf wax",
    "leaf-wax",
    "plant wax",
    "n-alkane",
    "n alkane",
    "n-alkanes",
    "alkane",
    "alkanoic",
    "fatty acid",
    "lipid",
    "biomarker",
    "biomarkers",
    "hydrogen isotope",
    "delta d",
    "δd",
    "d/h",
    "deuterium",
    "hydroclimate",
    "paleoclimate",
    "paleohydrolog",
    "soil",
    "sediment",
    "lacustrine",
    "lake",
    "loess",
    "terrestrial",
    "vegetation",
    "speleothem",
]

# Multiple query forms to try in case the first doesn't surface the right paper.
EXTRA_QUERY_FORMS = [
    "{author} {year} leaf wax n-alkane",
    "{author} {year} hydrogen isotope leaf wax",
    "{author} {year} n-alkane delta D precipitation",
    "{author} {year} plant wax dD",
    "{author} {year} alkane hydrogen isotope soil sediment",
    "{author} {year} leaf wax dD precipitation",
    "{author} {year} compound specific hydrogen isotope plant",
    "{author} {year} n-alkane lake sediment hydrogen",
]


def parse_author_year(label: str) -> tuple[str, str] | None:
    """Extract first author surname and year from a citation-like label.

    Handles:
      "Bai et al., 2011"
      "Polissar and Freeman 2010"
      "Polissar & Freeman, 2010"
      "Aichner et al 2019"
      "Ladd 2021 compilation (Herrmann et al., 2017)"
      "Ladd 2021 compilation (van der Veen et al., 2020)"
    """
    label = label.strip()
    m = re.search(r"\(([^)]+)\)", label)
    if m:
        label = m.group(1).strip()

    year_match = re.search(r"(19|20)\d{2}", label)
    if not year_match:
        return None
    year = year_match.group(0)

    # Strip the year and trailing 'in prep' style noise to find author chunk
    pre_year = label[: year_match.start()].strip().rstrip(",").strip()
    # Multi-word surname handling: "van der Veen et al" -> first 3 tokens before "et al"
    pre_year = re.sub(r"\bet\s+al\.?", "", pre_year, flags=re.IGNORECASE).strip()
    pre_year = pre_year.rstrip(",").strip()
    # If "and"/"&" present, take only the first surname
    pre_year = re.split(r"\s+(?:and|&)\s+", pre_year)[0].strip()
    if not pre_year:
        return None
    # Use the LAST token as the surname (handles "van der Veen" -> "Veen")
    tokens = pre_year.split()
    surname = tokens[-1] if tokens else None
    if not surname:
        return None
    return surname, year


def crossref_query(query: str, rows_n: int = 5) -> list[dict]:
    """Query Crossref bibliographic search."""
    params = {
        "query.bibliographic": query,
        "rows": str(rows_n),
    }
    url = f"{CROSSREF_URL}?{urllib.parse.urlencode(params)}"
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    try:
        with urllib.request.urlopen(req, timeout=30) as resp:
            data = json.loads(resp.read().decode("utf-8"))
    except Exception as exc:
        print(f"  Crossref error for {query!r}: {exc}", file=sys.stderr)
        return []
    return data.get("message", {}).get("items", [])


def crossref_author_filter_query(
    author: str, year: str, topic: str, rows_n: int = 10
) -> list[dict]:
    """Query Crossref with author filter and year filter (much more precise)."""
    params = {
        "query.author": author,
        "query.bibliographic": topic,
        "filter": f"from-pub-date:{year}-01-01,until-pub-date:{year}-12-31",
        "rows": str(rows_n),
    }
    url = f"{CROSSREF_URL}?{urllib.parse.urlencode(params)}"
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    try:
        with urllib.request.urlopen(req, timeout=30) as resp:
            data = json.loads(resp.read().decode("utf-8"))
    except Exception as exc:
        print(f"  Crossref filter error for author={author} year={year}: {exc}", file=sys.stderr)
        return []
    return data.get("message", {}).get("items", [])


REQUIRED_TOPIC_TERMS = [
    "leaf wax",
    "leaf-wax",
    "plant wax",
    "plant-wax",
    "n-alkane",
    "n alkane",
    "n-alkanes",
    "alkanoic",
    "alkanoate",
    "δd",
    "delta d",
    "deuterium",
    "hydrogen isotope",
    "d/h",
    "d2h",
    "δ2h",
]

JUNK_TITLES = ["editorial", "introduction", "preface", "errata", "erratum", "corrigendum"]


def score_item(item: dict, author: str, year: str) -> tuple[int, str]:
    """Score a Crossref item and return (score, reason).

    Strict rules:
    - Year must match exactly (no near-year acceptance — too risky for ambiguous names)
    - First author surname must match
    - Title must contain at least one REQUIRED_TOPIC_TERM (leaf wax / n-alkane / hydrogen isotope)
    - Title cannot be a junk word like "Editorial" or "Introduction"
    """
    reasons = []
    title = " ".join(item.get("title", []) or []).strip().lower()
    if not title:
        return -100, "no-title"
    if title in JUNK_TITLES:
        return -100, f"junk-title:{title}"

    # Year match — exact only
    issued = item.get("issued", {}).get("date-parts", [[None]])
    pub_year = None
    if issued and issued[0]:
        pub_year = issued[0][0]
    if not pub_year:
        return -100, "no-pub-year"
    if str(pub_year) != year:
        return -100, f"year-mismatch({pub_year}!={year})"
    reasons.append("year-exact")
    score = 5

    # Author surname must match first author (strict)
    authors = item.get("author", []) or []
    surnames = [(a.get("family") or "").lower() for a in authors]
    surnames = [s for s in surnames if s]
    a_lower = author.lower()
    if not surnames:
        return -100, "no-authors-on-record"
    if surnames[0] == a_lower:
        score += 8
        reasons.append("first-author-match")
    elif any(a_lower == s for s in surnames):
        # Acceptable but weaker
        score += 3
        reasons.append("author-in-list")
    elif any(a_lower in s.split() for s in surnames):
        # Substring match on a token (handles "Veen" inside "van der Veen" — though Crossref usually stores it as "Veen")
        score += 2
        reasons.append("author-substring")
    else:
        return -100, f"author-mismatch({a_lower} not in {surnames[:3]})"

    # REQUIRED topical relevance — must have at least one hard topic term
    required_hits = sum(1 for kw in REQUIRED_TOPIC_TERMS if kw in title)
    if required_hits == 0:
        return -100, f"no-required-topic title={title[:60]}"
    score += required_hits * 2
    reasons.append(f"required-topic-{required_hits}")

    # Bonus topical relevance
    topical_hits = sum(1 for kw in TOPICAL_KEYWORDS if kw in title)
    if topical_hits:
        score += topical_hits
        reasons.append(f"topic-{topical_hits}")

    return score, ",".join(reasons)


def find_doi(author: str, year: str) -> tuple[str | None, str | None, str | None]:
    """Try multiple Crossref queries; return (doi, title, reason). Returns
    (None, None, None) if no plausible match found."""
    seen_dois: set[str] = set()
    candidates: list[tuple[int, str, str, str]] = []  # (score, doi, title, reason)

    # Phase 1: bibliographic queries
    for form in EXTRA_QUERY_FORMS:
        query = form.format(author=author, year=year)
        items = crossref_query(query, rows_n=5)
        time.sleep(1.0)
        for item in items:
            doi = item.get("DOI")
            if not doi or doi in seen_dois:
                continue
            seen_dois.add(doi)
            score, reason = score_item(item, author, year)
            title = " ".join(item.get("title", []) or [])
            candidates.append((score, doi, title, reason))
        if candidates and max(c[0] for c in candidates) >= 14:
            break

    # Phase 2: if no strong match, try author+year filter searches with topic
    if not candidates or max(c[0] for c in candidates) < 14:
        for topic in [
            "leaf wax n-alkane hydrogen isotope",
            "n-alkane delta D plant",
            "leaf wax dD precipitation",
            "compound specific hydrogen lipid sediment",
        ]:
            items = crossref_author_filter_query(author, year, topic, rows_n=10)
            time.sleep(1.0)
            for item in items:
                doi = item.get("DOI")
                if not doi or doi in seen_dois:
                    continue
                seen_dois.add(doi)
                score, reason = score_item(item, author, year)
                title = " ".join(item.get("title", []) or [])
                candidates.append((score, doi, title, reason))
            if candidates and max(c[0] for c in candidates) >= 14:
                break

    if not candidates:
        return None, None, None
    candidates.sort(key=lambda c: c[0], reverse=True)
    best_score, best_doi, best_title, best_reason = candidates[0]
    # Strict threshold: at minimum year-exact (5) + first-author-match (8) + required-topic-1 (2) = 15
    # OR year-exact (5) + author-in-list (3) + required-topic-2 (4) + topic-2 (2) = 14
    if best_score < 14:
        return None, None, f"best-score={best_score} ({best_reason}) title={best_title[:70]}"
    return best_doi, best_title, f"score={best_score} {best_reason}"


# Sources for which there is no published DOI (in-prep / unpublished)
EXPECTED_BLANK = {
    "Ladd 2021 compilation (Ladd et al., SPCZ in prep)",
    "Ladd 2021 compilation (Ladd, Maloney et al (UW))",
}

# Post-Crossref corrections: each entry here was verified by direct Crossref
# search with location/context disambiguation above. These override whatever
# the automatic lookup found, or provide a DOI for a source the automatic
# lookup couldn't resolve. Do NOT use this table to fabricate DOIs — each
# entry is a Crossref-resolved DOI that the scoring missed or mis-matched.
MANUAL_VERIFIED_DOI: dict[str, str | None] = {
    # Mugler 2008 OG - Nam Co Tibet / Holzmaar (HZM in data)
    "Mugler et al 2008": "10.1016/j.orggeochem.2008.02.008",
    "Ladd 2021 compilation (Mugler et al 2008)": "10.1016/j.orggeochem.2008.02.008",
    # Freimuth 2019 GCA - Adirondack lake sediments (12 rows)
    "Ladd 2021 compilation (Freimuth et al 2019b)": "10.1016/j.gca.2019.08.026",
    # Freimuth 2019 OG - Brown's Lake Bog Ohio (1 row)
    "Ladd 2021 compilation (Freimuth et al., 2019a)": "10.1016/j.orggeochem.2019.01.006",
    # Peterse 2009 Biogeosciences - Soil n-alkane dD
    "Ladd 2021 compilation (Peterse et al., 2009)": "10.5194/bg-6-2799-2009",
    # Nelson 2018 - Ladd, Nelson, Schubert GCA - Lipid compound classes
    "Ladd 2021 compilation (Nelson et al., 2018)": "10.1016/j.gca.2018.06.005",
    # Wang 2017 heating/maturation paper (10.1016/j.orggeochem.2017.07.006) is
    # NOT the right paper for the 33 SW-Yunnan soil samples (22.2N, 100.5E).
    # No clear Crossref match found -> leave blank for manual resolution.
    "Wang et al., 2017": None,
    "Ladd 2021 compilation (Wang et al., 2017)": None,
    # Aichner 2019 Karakul Tajikistan: the 2018 CP match is a different
    # Aichner paper (Poland Younger Dryas). No clear 2019 Karakul DOI found
    # in Crossref. Leave blank for manual resolution.
    "Aichner et al 2019": None,
    # Huang 2016: Crossref best match (10.1016/j.quaint.2015.04.029) has
    # Li/Yang/Wang as first three authors — not Huang. No Huang-first 2016
    # leaf wax paper found in Crossref. Blank for manual resolution.
    "Huang et al., 2016": None,
    "Ladd 2021 compilation (Huang et al., 2016)": None,
    # Wang 2023 STOTEN paper (10.1016/j.scitotenv.2023.162970) has
    # Liu/Wang/Wang as first three authors — not Wang. Blank for manual.
    "Wang et al., 2023": None,
    # Seki et al 2010: the correct first-author paper is 10.1016/j.gca.2009.10.025
    # "A compound-specific n-alkane δ13C and δD approach for assessing source
    # and delivery processes of terrestrial organic matter". Crossref lists
    # the issued date as 2009 so strict year-exact rejected it; verified manually.
    "Seki et al., 2010": "10.1016/j.gca.2009.10.025",
}


def main() -> None:
    with CSV_PATH.open() as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        original_cols = reader.fieldnames

    print(f"Loaded {len(rows)} rows with columns: {original_cols}")

    # Build cache: source -> DOI from existing valid DOIs (excluding Ladd compilation)
    cache: dict[str, str] = {}
    for r in rows:
        src = r["source"]
        doi = r["DOI"].strip()
        if not doi:
            continue
        # Ladd's own paper keeps Ladd's DOI
        if src == LADD_PAPER_SOURCE:
            cache[src] = doi
            continue
        # Skip Ladd-compilation rows that have the Ladd DOI — those need replacement
        if doi == LADD_COMPILATION_DOI:
            continue
        if src not in cache:
            cache[src] = doi

    print(f"Sources already with DOI (kept as-is): {len(cache)}")

    all_sources = sorted({r["source"] for r in rows})
    needs_lookup = [
        s for s in all_sources if s not in cache and s not in EXPECTED_BLANK
    ]
    print(f"Sources requiring Crossref lookup: {len(needs_lookup)}")

    unresolved: list[tuple[str, str]] = []  # (source, reason)
    found_log: list[tuple[str, str, str]] = []  # (source, doi, title)

    for src in needs_lookup:
        parsed = parse_author_year(src)
        if not parsed:
            print(f"  SKIP (unparseable): {src!r}")
            unresolved.append((src, "unparseable"))
            continue
        author, year = parsed
        print(f"  Lookup: {author} {year}  ({src!r})")
        doi, title, reason = find_doi(author, year)
        if doi:
            print(f"    FOUND -> {doi}  [{reason}]")
            print(f"            title: {title[:90]}")
            cache[src] = doi
            found_log.append((src, doi, title))
        else:
            print(f"    NO MATCH  [{reason}]")
            unresolved.append((src, reason or "no-candidates"))

    # Apply manual verified overrides (may explicitly blank out a Crossref
    # mis-match, or supply a DOI the lookup missed).
    for src, verified in MANUAL_VERIFIED_DOI.items():
        if verified is None:
            cache.pop(src, None)
        else:
            cache[src] = verified

    # Apply DOIs to rows; drop columns
    out_cols = KEPT_COLS
    out_rows = []
    for r in rows:
        src = r["source"]
        new_doi = cache.get(src, "")
        # Special-case: any "Ladd compilation" row that still has the Ladd DOI but
        # is not the Ladd paper itself — replace with original ref DOI (or blank)
        if src != LADD_PAPER_SOURCE and r["DOI"].strip() == LADD_COMPILATION_DOI:
            r["DOI"] = new_doi
        else:
            r["DOI"] = new_doi if new_doi else r["DOI"].strip()
        out = {col: r.get(col, "") for col in out_cols}
        out_rows.append(out)

    with CSV_PATH.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=out_cols)
        writer.writeheader()
        writer.writerows(out_rows)

    print()
    print("=" * 60)
    print("VERIFICATION")
    print("=" * 60)
    print(f"Rows: {len(out_rows)}")
    unique_dois = {r["DOI"] for r in out_rows if r["DOI"]}
    print(f"Unique DOIs: {len(unique_dois)}")
    missing_doi_rows = [r for r in out_rows if not r["DOI"]]
    print(f"Rows with missing DOI: {len(missing_doi_rows)}")

    with CSV_PATH.open() as f:
        verify_reader = csv.DictReader(f)
        actual_cols = verify_reader.fieldnames
    print(f"Final columns: {actual_cols}")
    assert actual_cols == KEPT_COLS, f"Column mismatch: {actual_cols} != {KEPT_COLS}"
    print("Column list matches expected.")

    print()
    print("Top 10 unique sources by row count and assigned DOI:")
    from collections import Counter

    counts = Counter(r["source"] for r in out_rows)
    src_to_doi: dict[str, str] = {}
    for r in out_rows:
        src_to_doi.setdefault(r["source"], r["DOI"])
    for src, n in counts.most_common(10):
        doi = src_to_doi.get(src, "") or "(blank)"
        print(f"  [{n:4d}] {src!r}")
        print(f"         -> {doi}")

    print()
    print("All assigned DOIs (Crossref-found this run):")
    for src, doi, title in found_log:
        print(f"  {src!r}")
        print(f"    -> {doi}")
        print(f"       {title[:90]}")

    print()
    print(f"UNRESOLVED references ({len(unresolved)}):")
    for src, reason in unresolved:
        print(f"  - {src!r}  [reason: {reason}]")

    print()
    print(f"Sources expected to be blank (in-prep, unresolved by design):")
    for s in sorted(EXPECTED_BLANK):
        print(f"  - {s!r}")


if __name__ == "__main__":
    main()
