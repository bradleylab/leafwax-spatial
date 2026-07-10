import json, re
from pathlib import Path

data = json.loads(Path("manuscript/latex/_parsed.json").read_text())
parsed = data["parsed"]

# Cleanup but preserve trailing period on journal abbreviations
for p in parsed:
    for f in ("publisher", "address", "title"):
        if f in p and p[f]:
            v = p[f].strip()
            v = re.sub(r"\*", "", v)
            v = v.strip().rstrip(",.")
            p[f] = v
    # journal: strip * and trailing comma but KEEP trailing period
    if "journal" in p and p["journal"]:
        v = p["journal"].strip()
        v = re.sub(r"\*", "", v).strip().rstrip(",")
        # If it looks like an abbreviation (≥1 dot inside) ensure terminal period
        if re.search(r"\.", v) and not v.endswith("."):
            v = v + "."
        p["journal"] = v

# Special-case overrides
override = {
    "friedl2019": {"publisher": "NASA EOSDIS Land Processes Distributed Active Archive Center"},
    "wood2017": {"address": "Boca Raton, FL"},
    "box2013": {"address": "Hoboken, NJ"},
    "stein1999": {"address": "New York, NY"},
    "diggle2007": {"address": "New York, NY"},
    "bolker2008": {"address": "Princeton, NJ"},
    "jaccard2003": {"address": "Thousand Oaks, CA"},
}
for k, fields in override.items():
    for p in parsed:
        if p["key"] == k:
            for kk, vv in fields.items():
                p[kk] = vv

new_cites = [
    {"key": "khan2026", "type": "article", "author": "Khan, K. and Berrett, M.", "year": "2026",
     "title": "Re-Thinking Spatial Confounding in Spatial Linear Mixed Models",
     "journal": "Statistical Science", "volume": "41", "number": "2", "doi": "10.1214/25-sts976"},
    {"key": "khan2022", "type": "article", "author": "Khan, K. and Calder, C. A.", "year": "2022",
     "title": "Restricted Spatial Regression Methods: Implications for Inference",
     "journal": "Journal of the American Statistical Association",
     "volume": "117", "number": "537", "pages": "482--495", "doi": "10.1080/01621459.2020.1788949"},
    {"key": "gilbert2024", "type": "article", "author": "Gilbert, B. and Datta, A. and Casey, J. A. and Ogburn, E. L.", "year": "2024",
     "title": "Consistency of common spatial estimators under spatial confounding",
     "journal": "Biometrics", "volume": "81", "number": "1"},
    {"key": "baan2025", "type": "article", "author": "Baan, J. and Holloway-Phillips, M. and Nelson, D. B. and de Vos, R. C. H. and Kahmen, A.", "year": "2025",
     "title": "Phylogenetic and biochemical drivers of plant species variation in organic compound hydrogen stable isotopes",
     "journal": "New Phytologist", "volume": "246", "doi": "10.1111/nph.20430"},
    {"key": "polissar2025", "type": "article", "author": "Polissar, P. J. and Karp, A. T. and D'Andrea, W. J.", "year": "2025",
     "title": "Mixed messages: Unmixing sedimentary molecular distributions reveals source contributions and isotopic values",
     "journal": "Geochimica et Cosmochimica Acta", "volume": "380", "doi": "10.1016/j.gca.2025.03.001"},
    {"key": "sivula2020", "type": "article", "author": "Sivula, T. and Magnusson, M. and Vehtari, A.", "year": "2020",
     "title": "Uncertainty in Bayesian Leave-One-Out Cross-Validation Based Model Comparison",
     "journal": "arXiv preprint", "note": "arXiv:2008.10296"},
]
existing_keys = {p["key"] for p in parsed}
for c in new_cites:
    if c["key"] not in existing_keys:
        parsed.append(c)

def escape_title(t):
    if not t: return t
    repl = {
        "δ": r"$\delta$", "Δ": r"$\Delta$", "ε": r"$\varepsilon$",
        "α": r"$\alpha$", "β": r"$\beta$", "λ": r"$\lambda$",
        "ρ": r"$\rho$", "σ": r"$\sigma$", "μ": r"$\mu$",
        "²": r"$^{2}$", "³": r"$^{3}$",
        "₂": r"$_{2}$", "₃": r"$_{3}$", "₄": r"$_{4}$",
        "−": "-", "–": "--", "—": "---",
        "’": "'", "‘": "'", "“": "``", "”": "''",
    }
    for k, v in repl.items():
        t = t.replace(k, v)
    # Brace-protect ALL-CAPS+ tokens (≥2 caps) for case preservation under elsarticle-harv
    t = re.sub(r"\b[A-Z]{2,}[A-Za-z0-9]*\b", lambda m: "{" + m.group(0) + "}", t)
    return t

def escape_field(v):
    if not v: return v
    v = str(v)
    v = v.replace("–", "--").replace("—", "---")
    v = v.replace("’", "'").replace("‘", "'")
    # Escape & for BibTeX (publisher: "Chapman & Hall/CRC" -> "Chapman \& Hall/CRC")
    v = v.replace("&", r"\&")
    return v

def emit(p):
    t = p["type"]
    fields = {}
    for f in ("author", "title", "journal", "year", "volume", "number", "pages", "publisher", "address", "edition", "doi", "note"):
        if f in p and p[f]:
            if f == "title":
                fields[f] = escape_title(p[f])
            elif f == "author":
                a = str(p[f]).replace("’", "'").replace("‘", "'").replace("–", "-")
                fields[f] = a
            else:
                fields[f] = escape_field(p[f])
    lines = [f"@{t}{{{p['key']},"]
    for f, v in fields.items():
        lines.append(f"  {f} = {{{v}}},")
    lines[-1] = lines[-1].rstrip(",")
    lines.append("}")
    return "\n".join(lines)

bib_text = "% Auto-generated from manuscript/drafts/leafwax_gca_text_2026-05-05.md and supplement\n"
bib_text += "% Generated 2026-05-05; verify DOIs and missing fields against Zotero before submission\n\n"
parsed_sorted = sorted(parsed, key=lambda x: x["key"])
bib_text += "\n\n".join(emit(p) for p in parsed_sorted)
bib_text += "\n"

Path("manuscript/latex/references.bib").write_text(bib_text)
print(f"wrote {len(parsed_sorted)} entries to manuscript/latex/references.bib")

# Save final keys
keys = [p["key"] for p in parsed_sorted]
Path("manuscript/latex/_keys.json").write_text(json.dumps(keys, indent=2))
