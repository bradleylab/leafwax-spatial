import re, json, unicodedata
from pathlib import Path

main_md = Path("manuscript/drafts/leafwax_gca_text_2026-05-05.md").read_text()
supp_md = Path("manuscript/drafts/leafwax_gca_supplement_2026-05-05.md").read_text()

def split_entries(block):
    return [c.strip() for c in re.split(r"\n\s*\n", block.strip()) if c.strip() and re.match(r"^[A-ZÀ-ÿ]", c.strip())]

main_refs = split_entries(main_md.split("**References**", 1)[1])
supp_refs = split_entries(supp_md.split("**Supplementary References**", 1)[1])

def primary_key(entry):
    yr_m = re.search(r"\((\d{4})([a-z]?)\)", entry)
    if not yr_m: return None
    year, suffix = yr_m.group(1), yr_m.group(2) or ""
    sur_m = re.match(r"^([A-ZÀ-Ö][A-Za-zÀ-ÿ\-']+)", entry)
    if not sur_m: return None
    sur = unicodedata.normalize("NFKD", sur_m.group(1)).encode("ascii", "ignore").decode("ascii").lower()
    sur = re.sub(r"[^a-z]", "", sur)
    return f"{sur}{year}{suffix}"

dedup = {}
for e in main_refs + supp_refs:
    k = primary_key(e)
    if k and k not in dedup:
        dedup[k] = e

INITIAL_RE = re.compile(r"^[A-ZÀ-Ö](?:\.\-?[A-ZÀ-Ö])?\.?$")
INITIAL_TOKEN_RE = re.compile(r"^[A-ZÀ-Ö](?:\.[A-ZÀ-Ö])*\.?$")
HYPHEN_INITIAL_RE = re.compile(r"^[A-ZÀ-Ö]\.?\-[A-ZÀ-Ö]\.?$")

def normalize_initials(seq):
    """Take a sequence of tokens that are all initials; return canonical form like 'X. Y.' or 'X.-Y.'"""
    out = []
    for t in seq:
        # Hyphenated initials: 'X-Y' or 'X.-Y.' or 'X.-Y'
        if "-" in t:
            parts = t.split("-")
            normed = "-".join(p.rstrip(".") + "." for p in parts if p)
            out.append(normed)
        else:
            # Strip dots, then re-add per letter
            letters = re.sub(r"[\.\s]", "", t)
            if letters:
                out.append(" ".join(c + "." for c in letters))
    return " ".join(out)

def parse_authors(pre_year):
    s = pre_year.strip().rstrip(",.")
    # 'and' -> ';'
    s = re.sub(r"\s+and\s+", " ; ", s)
    # split by commas
    raw_segments = [p.strip() for p in s.split(",") if p.strip()]
    # combine continuations: a segment that starts with an initial token belongs to the previous author
    authors_combined, cur = [], []
    for p in raw_segments:
        toks = p.split()
        if cur and toks and (INITIAL_RE.match(toks[0]) or HYPHEN_INITIAL_RE.match(toks[0])):
            cur.append(p)
        else:
            if cur: authors_combined.append(", ".join(cur))
            cur = [p]
    if cur: authors_combined.append(", ".join(cur))
    expanded = []
    for a in authors_combined:
        for sub in a.split(" ; "):
            sub = sub.strip()
            if sub: expanded.append(sub)
    bib = []
    for a in expanded:
        toks = a.split()
        if not toks: continue
        # find boundary: first initial-looking token
        boundary = None
        for i, t in enumerate(toks):
            if INITIAL_RE.match(t) or HYPHEN_INITIAL_RE.match(t) or INITIAL_TOKEN_RE.match(t):
                boundary = i
                break
        if boundary is None:
            bib.append(a)
        elif boundary == 0:
            # the entire thing is initials — odd, just emit raw
            bib.append(a)
        else:
            surname = " ".join(toks[:boundary])
            initials = normalize_initials(toks[boundary:])
            bib.append(f"{surname}, {initials}" if initials else surname)
    return " and ".join(bib)

def parse_entry(text):
    text = text.strip()
    yr_m = re.search(r"\((\d{4})([a-z]?)\)", text)
    if not yr_m: return None
    year = yr_m.group(1)
    pre_year = text[:yr_m.start()].strip().rstrip(",.")
    post_year = text[yr_m.end():].strip()
    bib_authors = parse_authors(pre_year)
    post_year_stripped = post_year.lstrip()
    title_italic = post_year_stripped.startswith("*")

    if title_italic:
        m = re.search(r"\*([^*]+)\*", post_year)
        if not m:
            return {"key": primary_key(text), "type": "misc", "author": bib_authors, "year": year, "title": post_year.rstrip(".")}
        title = m.group(1).strip().rstrip(".,")
        after = post_year[m.end():].strip()
        ed_m = re.search(r"(\d+(?:st|nd|rd|th))\s*ed\.?", after)
        edition = ed_m.group(1) if ed_m else None
        # strip leading punctuation/whitespace and edition
        after_clean = re.sub(r"^[\.,\s]+", "", after)
        after_clean = re.sub(r"^\d+(?:st|nd|rd|th)\s*ed\.?[,\.]?\s*", "", after_clean).strip().rstrip(".")
        chunks = [c.strip() for c in after_clean.split(",") if c.strip() and c.strip() != "."]
        publisher = chunks[0] if chunks else None
        address = ", ".join(chunks[1:]) if len(chunks) > 1 else None
        out = {"key": primary_key(text), "type": "book", "author": bib_authors, "year": year, "title": title}
        if publisher: out["publisher"] = publisher
        if address: out["address"] = address
        if edition: out["edition"] = edition
        return out

    m = re.search(r"\*([^*]+)\*", post_year)
    if m:
        title = post_year[:m.start()].strip().rstrip(".,")
        journal = m.group(1).strip().rstrip(",")
        rest = post_year[m.end():].strip()
    else:
        return {"key": primary_key(text), "type": "misc", "author": bib_authors, "year": year, "title": post_year.rstrip(".")}

    vol_m = re.search(r"\*\*(\d+(?:[-–]\d+)?)\*\*", rest)
    volume = vol_m.group(1).replace("–", "--") if vol_m else None
    page_m = re.search(r",\s*([\dee\w\.\s\-–]+?)\.?\s*$", rest)
    pages = None
    if page_m:
        p = page_m.group(1).strip().rstrip(".").replace("–", "--")
        if re.match(r"^([\dee]+([\-–\.]\d+)*|[a-z]?\d+[a-z]?(\.\d+)?)", p):
            pages = p
    fields = {"key": primary_key(text), "type": "article", "author": bib_authors, "year": year, "title": title}
    if journal: fields["journal"] = journal
    if volume: fields["volume"] = volume
    if pages: fields["pages"] = pages
    return fields

parsed = []
for k, e in dedup.items():
    p = parse_entry(e)
    if p: parsed.append(p)

from collections import Counter
print(Counter(p["type"] for p in parsed))

# check the previously broken ones
for k in ["vehtari2021", "wang2023", "national2006", "iaea2015", "bowen2018", "amatulli2018", "friedl2019", "danielson2011", "konecky2019a", "obrien2007", "feakins2013", "polissar2014", "cernusak2016"]:
    for p in parsed:
        if p["key"] == k:
            print(json.dumps(p, indent=2, ensure_ascii=False))

Path("manuscript/latex/_parsed.json").write_text(json.dumps({"parsed": parsed, "raw": dedup}, indent=2, ensure_ascii=False))
