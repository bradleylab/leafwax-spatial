"""Convert leafwax_gca_supplement_2026-05-05.md → supplement.tex (elsarticle, manual S-numbering)."""
import re, json, sys, unicodedata
from pathlib import Path

KEYS = set(json.loads(Path("manuscript/latex/_keys.json").read_text()))
SURNAME_ALIASES = {"council": "national", "wmo": "iaea"}
SURNAME_BLACKLIST = {"jr", "sr", "iii", "ii", "iv"}
UNMATCHED = []

def normalize_surname(s):
    s = unicodedata.normalize("NFKD", s).encode("ascii", "ignore").decode("ascii").lower()
    return re.sub(r"[^a-z]", "", s)

def lookup_key(surname, year_letter):
    n = normalize_surname(surname)
    if n in SURNAME_BLACKLIST: return None
    if n in SURNAME_ALIASES: n = SURNAME_ALIASES[n]
    cand = f"{n}{year_letter}"
    if cand in KEYS: return cand
    if re.match(r"\d{4}[a-z]$", year_letter):
        cand2 = f"{n}{year_letter[:-1]}"
        if cand2 in KEYS: return cand2
    return None

def parse_one_cite(s):
    s = s.strip().rstrip(",.").strip()
    s = re.sub(r"^(?:e\.g\.,?|see|cf\.,?|after|following|sensu|including)\s+", "", s, flags=re.IGNORECASE).strip()
    multi_m = re.match(r"^([A-ZÀ-Ö][A-Za-zÀ-ÿ\-'\s]+?(?:\s+et\s+al\.?)?)\s*,?\s*(\d{4}[a-z]?(?:\s*,\s*\d{4}[a-z]?)+)\s*$", s)
    if multi_m:
        prefix = multi_m.group(1).strip().rstrip(",")
        years = [y.strip() for y in re.split(r"\s*,\s*", multi_m.group(2))]
        first_word = re.match(r"^([A-ZÀ-Ö][A-Za-zÀ-ÿ\-']+)", prefix).group(1)
        keys = []
        for y in years:
            k = lookup_key(first_word, y)
            if k: keys.append(k)
            else: UNMATCHED.append(f"{prefix}, {y}"); return None
        return keys
    single_m = re.match(r"^(.+?),\s*(\d{4}[a-z]?)\s*$", s)
    if single_m:
        author_block = single_m.group(1).strip(); year = single_m.group(2)
        first = re.match(r"^([A-ZÀ-Ö][A-Za-zÀ-ÿ\-']+)", author_block)
        if not first: return None
        k = lookup_key(first.group(1), year)
        if k: return [k]
        UNMATCHED.append(f"{author_block}, {year}"); return None
    nc_m = re.match(r"^([A-ZÀ-Ö][A-Za-zÀ-ÿ\-']+(?:\s+(?:and|&)\s+[A-ZÀ-Ö][A-Za-zÀ-ÿ\-']+)?(?:\s+et\s+al\.?)?)\s+(\d{4}[a-z]?)\s*$", s)
    if nc_m:
        first = re.match(r"^([A-ZÀ-Ö][A-Za-zÀ-ÿ\-']+)", nc_m.group(1)).group(1)
        k = lookup_key(first, nc_m.group(2))
        if k: return [k]
        UNMATCHED.append(s)
    return None

def replace_paren_cites(text):
    pat = re.compile(r"\(([^()]*?\d{4}[a-z]?(?:[^()]*?)?)\)")
    out, i = [], 0
    while True:
        m = pat.search(text, i)
        if not m:
            out.append(text[i:]); break
        out.append(text[i:m.start()])
        body = m.group(1)
        if not re.search(r"\b\d{4}[a-z]?\b", body):
            out.append(m.group(0)); i = m.end(); continue
        parts = [p.strip() for p in body.split(";")]
        all_keys, prose_pieces, any_failure = [], [], False
        for p in parts:
            k = parse_one_cite(p)
            if k: all_keys.extend(k)
            else: prose_pieces.append(p); any_failure = True
        if all_keys and not any_failure:
            first = parts[0]
            pm = re.match(r"^(e\.g\.,?|see|cf\.,?|after|following|sensu)\s+", first, flags=re.IGNORECASE)
            if pm:
                pre = pm.group(1).rstrip(",.") + ","
                out.append(f"\\citep[{pre}][]{{{','.join(all_keys)}}}")
            else:
                out.append(f"\\citep{{{','.join(all_keys)}}}")
        elif all_keys and prose_pieces:
            rebuilt = []
            for p in parts:
                k = parse_one_cite(p)
                if k: rebuilt.append(f"\\citealp{{{','.join(k)}}}")
                else: rebuilt.append(p)
            out.append("(" + "; ".join(rebuilt) + ")")
        else:
            out.append(m.group(0))
        i = m.end()
    return "".join(out)

def replace_narrative_cites(text):
    def repl(m):
        block = m.group(1); year = m.group(2)
        first = re.match(r"^([A-ZÀ-Ö][A-Za-zÀ-ÿ\-']+)", block).group(1)
        k = lookup_key(first, year)
        if k: return f"\\citet{{{k}}}"
        UNMATCHED.append(f"{block} ({year})")
        return m.group(0)
    pat = re.compile(r"([A-ZÀ-Ö][A-Za-zÀ-ÿ\-']+(?:\s+(?:and|&)\s+[A-ZÀ-Ö][A-Za-zÀ-ÿ\-']+)?(?:\s+et\s+al\.?)?)\s*\((\d{4}[a-z]?)\)")
    return pat.sub(repl, text)

MATH_SENTINEL = "MMATHPROTECTNN"  # printable, unlikely to occur in prose

def md_inline(s):
    # Protect inline math $...$ from markdown processing
    math_blocks = []
    def stash_math(m):
        math_blocks.append(m.group(0))
        return f"{MATH_SENTINEL}{len(math_blocks)-1}{MATH_SENTINEL}"
    s = re.sub(r"\$[^$\n]+?\$", stash_math, s)

    s = s.replace(r"\'", "'").replace(r"\<", "<").replace(r"\>", ">").replace(r"\|", "|").replace(r"\~", "~").replace(r"\*", "*")
    s = s.replace("ε~app~", r"$\varepsilon_{\mathrm{app}}$")
    s = s.replace("²ε~app~", r"$^{2}\varepsilon_{\mathrm{app}}$")
    s = re.sub(r"([A-Za-z])̂", r"$\\hat{\1}$", s)
    s = re.sub(r"\^([^\s^]+?)\^", r"\\textsuperscript{\1}", s)
    s = re.sub(r"~([^\s~]+?)~", r"\\textsubscript{\1}", s)
    s = re.sub(r"\*\*(.+?)\*\*", r"\\textbf{\1}", s)
    s = re.sub(r"(?<!\\)\*([^*\n]+?)\*", r"\\emph{\1}", s)
    repl = {
        "δ": r"$\delta$", "Δ": r"$\Delta$", "ε": r"$\varepsilon$",
        "α": r"$\alpha$", "β": r"$\beta$", "λ": r"$\lambda$",
        "ρ": r"$\rho$", "σ": r"$\sigma$", "μ": r"$\mu$", "τ": r"$\tau$",
        "χ": r"$\chi$", "ν": r"$\nu$", "η": r"$\eta$", "θ": r"$\theta$",
        "ϕ": r"$\phi$", "φ": r"$\varphi$", "γ": r"$\gamma$",
        "‰": r"\textperthousand{}",
        "⁰": r"$^{0}$", "¹": r"$^{1}$", "²": r"$^{2}$", "³": r"$^{3}$",
        "⁴": r"$^{4}$", "⁵": r"$^{5}$", "⁶": r"$^{6}$", "⁷": r"$^{7}$",
        "⁸": r"$^{8}$", "⁹": r"$^{9}$",
        "₀": r"$_{0}$", "₁": r"$_{1}$", "₂": r"$_{2}$", "₃": r"$_{3}$",
        "₄": r"$_{4}$", "₅": r"$_{5}$", "₆": r"$_{6}$", "₇": r"$_{7}$",
        "₈": r"$_{8}$", "₉": r"$_{9}$",
        "ᵢ": r"$_{i}$", "ⱼ": r"$_{j}$", "ₖ": r"$_{k}$", "ₙ": r"$_{n}$",
        "−": "-", "–": "--", "—": "---", "‐": "-",
        "≈": r"$\approx$", "≤": r"$\leq$", "≥": r"$\geq$",
        "×": r"$\times$", "·": r"$\cdot$",
        "⁻": r"$^{-}$", "⁺": r"$^{+}$",
        "°": r"\textdegree{}",
        "→": r"$\rightarrow$", "←": r"$\leftarrow$",
        "’": "'", "‘": "'", "“": "``", "”": "''",
        "…": r"\ldots{}",
        "Å": r"\AA{}",
        "%": r"\%",
        "±": r"$\pm$",
        "√": r"$\sqrt{}$",
        "∑": r"$\sum$", "∫": r"$\int$", "∞": r"$\infty$",
        "∝": r"$\propto$", "∈": r"$\in$",
    }
    for k, v in repl.items():
        s = s.replace(k, v)
    s = re.sub(r"(?<!\$)~(?=[\d])", r"$\\sim$", s)
    s = re.sub(r"(?<!\$)~(?=[A-Za-z])", r"$\\sim$", s)
    s = re.sub(r"<([A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,})>", r"\\href{mailto:\1}{\1}", s)
    s = re.sub(r"(?<![\\$])<", r"$<$", s)
    s = re.sub(r"(?<![\\$])>", r"$>$", s)

    # Escape _, #, & outside math (but math is currently stashed as sentinels — so safe everywhere)
    out = []
    i = 0
    while i < len(s):
        c = s[i]
        if c == '$':  # math toggling — but actual $ ... $ are stashed; only stray $s remain
            out.append(c); i += 1; continue
        if c == '_':
            if i > 0 and s[i-1] == '\\':
                out.append(c); i += 1; continue
            out.append(r'\_'); i += 1; continue
        if c == '#':
            out.append(r'\#'); i += 1; continue
        if c == '&':
            if i > 0 and s[i-1] == '\\':
                out.append(c); i += 1; continue
            out.append(r'\&'); i += 1; continue
        out.append(c); i += 1
    result = "".join(out)
    # Restore math
    for idx, mb in enumerate(math_blocks):
        result = result.replace(f"{MATH_SENTINEL}{idx}{MATH_SENTINEL}", mb)
    return result

def convert_text(text):
    text = replace_narrative_cites(text)
    text = replace_paren_cites(text)
    text = md_inline(text)
    return text

SRC = Path("manuscript/drafts/leafwax_gca_supplement_2026-05-05.md").read_text()
body_text = SRC.split("**Supplementary References**", 1)[0]

def scan_body(body):
    lines = body.splitlines()
    out, i, in_eq, eq_buf = [], 0, False, []
    HEADER = re.compile(r"^\*\*(S\d+(?:\.\d+)*)\.?\s+([^*]+?)\*\*\s*$")
    SPECIAL = re.compile(r"^\*\*([A-Z][A-Za-z][A-Za-z\s:]*?):?\*\*\s*$")
    while i < len(lines):
        line = lines[i]
        if line.strip().startswith("$$"):
            stripped = line.strip()
            if in_eq:
                tail = stripped[2:].strip().rstrip("$").strip()
                if tail: eq_buf.append(tail)
                out.append(("equation", "\n".join(s for s in eq_buf if s)))
                in_eq, eq_buf = False, []
            else:
                if stripped.endswith("$$") and len(stripped) > 4:
                    out.append(("equation", stripped[2:-2].strip()))
                else:
                    in_eq = True
                    inner = stripped[2:].strip()
                    if inner: eq_buf.append(inner)
            i += 1; continue
        if in_eq:
            eq_buf.append(line); i += 1; continue
        m = HEADER.match(line.strip())
        if m:
            num, title = m.group(1), m.group(2).strip().rstrip(":")
            depth = num.count(".") + 1
            out.append(("header", depth, title, num))
            i += 1; continue
        eql_m = re.match(r"^\(\s*Suppl(?:emental|ementary)\s+Equations?\s+(\d+)\s*\)\s*$", line.strip(), flags=re.IGNORECASE)
        if eql_m:
            out.append(("equation_label", eql_m.group(1))); i += 1; continue
        sp = SPECIAL.match(line.strip())
        if sp:
            out.append(("special", sp.group(1).strip().rstrip(":"))); i += 1; continue
        if re.match(r"^\s*[-*]\s+", line):
            items = []
            while i < len(lines) and re.match(r"^\s*[-*]\s+", lines[i]):
                items.append(re.sub(r"^\s*[-*]\s+", "", lines[i])); i += 1
            out.append(("list", items))
            continue
        out.append(("para", line))
        i += 1
    return out

elements = scan_body(body_text)
grouped, buf = [], []
for el in elements:
    if el[0] == "para":
        if el[1].strip() == "":
            if buf: grouped.append(("para", "\n".join(buf).strip())); buf = []
        else:
            buf.append(el[1])
    else:
        if buf: grouped.append(("para", "\n".join(buf).strip())); buf = []
        grouped.append(el)
if buf: grouped.append(("para", "\n".join(buf).strip()))

first_header_idx = None
for k, el in enumerate(grouped):
    if el[0] == "header" and first_header_idx is None:
        first_header_idx = k

body_chunks = []
for el in grouped[first_header_idx:]:
    t = el[0]
    if t == "header":
        depth, title, num = el[1], el[2], el[3]
        cmd_map = {1: r"\section*", 2: r"\subsection*", 3: r"\subsubsection*", 4: r"\paragraph"}
        cmd = cmd_map.get(depth, r"\paragraph")
        body_chunks.append(f"\n{cmd}{{{num}~~{md_inline(title)}}}\n\\addcontentsline{{toc}}{{section}}{{{num} {title}}}\n")
    elif t == "special":
        nm = el[1].lower()
        if nm in ("supplementary references",):
            break
        else:
            body_chunks.append(f"\n\\paragraph{{{el[1]}}}\n")
    elif t == "para":
        if el[1].strip().startswith("**Spatial modeling improves") or el[1].strip() == "**Supplementary Material**":
            continue
        body_chunks.append("\n" + convert_text(el[1]) + "\n")
    elif t == "equation":
        body_chunks.append("\n\\begin{equation}\n" + el[1] + "\n\\end{equation}\n")
    elif t == "equation_label":
        pass
    elif t == "list":
        body_chunks.append("\n\\begin{itemize}\n")
        for item in el[1]:
            body_chunks.append("  \\item " + convert_text(item) + "\n")
        body_chunks.append("\\end{itemize}\n")

body_tex = "".join(body_chunks)

fig_caps = {}
for m in re.finditer(r"\*\*Figure (S\d+)\.\s+(.*?)(?=\n\n\*\*Figure S\d+\.|\n\n\*\*Table |\n\n\*\*Supplementary References\*\*|\Z)", SRC, re.DOTALL):
    label = m.group(1); rest = m.group(2)
    mm = re.search(r"\.\*\*\s", rest)
    if mm:
        short = rest[:mm.end()-3].strip().rstrip("*").strip()
        full = rest[mm.end():].strip()
    else:
        short = rest[:80].strip(); full = rest[80:].strip()
    fig_caps[label] = (short, full)

tab_caps = {}
for m in re.finditer(r"\*\*Table (S\d+)\.\*\*\s+(.*?)(?=\n\n\*\*Table S\d+\.|\n\n\*\*Figure S\d+\.|\n\n\*\*Supplementary References\*\*|\Z)", SRC, re.DOTALL):
    tab_caps[m.group(1)] = m.group(2).strip()

preamble = r"""\documentclass[3p,times,authoryear]{elsarticle}

\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage{amsmath,amssymb}
\usepackage{graphicx}
\usepackage{booktabs}
\usepackage{longtable}
\usepackage{array}
\usepackage{xcolor,colortbl}
\usepackage{textcomp}
\usepackage{lscape}
\usepackage{lineno}
\usepackage{natbib}
\usepackage{url}
\usepackage{hyperref}

\modulolinenumbers[1]
\linenumbers

\graphicspath{{../figures/supplement_figs/}{../figures/main_figs/}}

\journal{Geochimica et Cosmochimica Acta}

\renewcommand{\thefigure}{S\arabic{figure}}
\renewcommand{\thetable}{S\arabic{table}}
\renewcommand{\theequation}{S\arabic{equation}}
"""

front = r"""
\begin{document}

\begin{frontmatter}
\title{Supplementary Material for ``Spatial modeling improves the calibration of leaf wax hydrogen isotopes to precipitation''}
\author[wustl]{Alexander S. Bradley}
\address[wustl]{Department of Earth, Environmental, and Planetary Sciences,
                Washington University in St. Louis,
                1 Brookings Drive, Saint Louis, Missouri 63130, USA}
\end{frontmatter}

"""

postamble = r"""

\bibliographystyle{elsarticle-harv}
\bibliography{references}

\end{document}
"""

fig_files = {
    "S1": "Figure_S1_spatial_weighting.pdf",
    "S2": "Figure_S2_pairwise_correlations.pdf",
    "S3": "Figure_S3_residuals.pdf",
    "S4": "Figure_S4_ols_spatial_scale.pdf",
    "S5": "Figure_S5_elevation.pdf",
}

fig_tex = ["\n\\clearpage\n"]
for n in sorted(fig_caps.keys()):
    short, full = fig_caps[n]
    fname = fig_files.get(n)
    if not fname: continue
    short_clean = convert_text(short)
    full_clean = convert_text(full)
    num = int(n[1:])
    fig_tex.append(rf"""
\setcounter{{figure}}{{{num-1}}}
\begin{{figure}}[p]
  \centering
  \includegraphics[width=0.95\linewidth]{{{fname}}}
  \caption[{short_clean}]{{\textbf{{{short_clean}.}} {full_clean}}}
  \label{{fig:fig{n.lower()}}}
\end{{figure}}
""")

tab_tex = ["\n\\clearpage\n"]

s1_full = Path("manuscript/tables/Table_S1_priors.tex").read_text()
s1_match = re.search(r"\\begin\{longtable\}.*?\\end\{longtable\}", s1_full, re.DOTALL)
if s1_match:
    s1_body = s1_match.group(0)
    tab_tex.append(f"""
\\setcounter{{table}}{{0}}
\\begin{{landscape}}
{s1_body}
\\end{{landscape}}
""")

s2_caption = tab_caps.get("S2", "Regional in-sample performance metrics for spatial models")
_s2_body_clean_path = Path("manuscript/latex/_table_s2_body_clean.tex")
_s2_caption_tex = convert_text(s2_caption)
_s2_block = (
    "\n\\clearpage\n"
    "\\setcounter{table}{1}\n"
    "\\begin{table}[p]\n"
    "  \\centering\n"
    "  \\caption{\\textbf{Table S2.} " + _s2_caption_tex + "}\n"
    "  \\label{tab:tabS2}\n"
    "  \\footnotesize\n"
    "  \\begin{tabular}{@{}llrrr@{}}\n"
    "    \\toprule\n"
    "    Model & Region & $n$ & RMSE (\\textperthousand) & $R^2$ \\\\\n"
    "    \\midrule\n"
    + _s2_body_clean_path.read_text() + "\n"
    "    \\bottomrule\n"
    "  \\end{tabular}\n"
    "\\end{table}\n"
)
tab_tex.append(_s2_block)
_DEAD_BLOCK_BELOW = """
foo bar baz placeholder
\rowcolor{lightgray} \textbf{baseline} & \textbf{Overall} & \textbf{1129} & \textbf{21.2} & \textbf{0.696} \\
 & Africa & 143 & 20.5 & 0.021 \\
 & Americas & 373 & 19.8 & 0.835 \\
 & Asia & 551 & 20.5 & 0.460 \\
 & Europe & 30 & 17.3 & 0.718 \\
 & Oceania & 32 & 44.1 & 0.685 \\
\rowcolor{lightgray} \textbf{baseline\_sp} & \textbf{Overall} & \textbf{1129} & \textbf{15.9} & \textbf{0.829} \\
 & Africa & 143 & 12.2 & 0.640 \\
 & Americas & 373 & 14.7 & 0.900 \\
 & Asia & 551 & 17.7 & 0.576 \\
 & Europe & 30 & 12.4 & 0.755 \\
 & Oceania & 32 & 15.1 & 0.828 \\
\rowcolor{lightgray} \textbf{baseline\_env} & \textbf{Overall} & \textbf{1129} & \textbf{20.3} & \textbf{0.724} \\
 & Africa & 143 & 19.4 & 0.119 \\
 & Americas & 373 & 19.4 & 0.845 \\
 & Asia & 551 & 20.0 & 0.470 \\
 & Europe & 30 & 14.3 & 0.741 \\
 & Oceania & 32 & 36.7 & 0.690 \\
\rowcolor{lightgray} \textbf{baseline\_env\_sp} & \textbf{Overall} & \textbf{1129} & \textbf{15.6} & \textbf{0.836} \\
 & Africa & 143 & 11.5 & 0.682 \\
 & Americas & 373 & 14.4 & 0.905 \\
 & Asia & 551 & 17.4 & 0.591 \\
 & Europe & 30 & 12.3 & 0.760 \\
 & Oceania & 32 & 15.0 & 0.825 \\
\rowcolor{lightgray} \textbf{baseline\_veg} & \textbf{Overall} & \textbf{1129} & \textbf{20.1} & \textbf{0.729} \\
 & Africa & 143 & 14.9 & 0.551 \\
 & Americas & 373 & 19.4 & 0.848 \\
 & Asia & 551 & 19.9 & 0.476 \\
 & Europe & 30 & 16.3 & 0.690 \\
 & Oceania & 32 & 41.7 & 0.661 \\
\rowcolor{lightgray} \textbf{baseline\_veg\_sp} & \textbf{Overall} & \textbf{1129} & \textbf{15.7} & \textbf{0.834} \\
 & Africa & 143 & 11.1 & 0.698 \\
 & Americas & 373 & 14.2 & 0.906 \\
 & Asia & 551 & 17.8 & 0.573 \\
 & Europe & 30 & 12.4 & 0.757 \\
 & Oceania & 32 & 14.7 & 0.836 \\
\rowcolor{lightgray} \textbf{full} & \textbf{Overall} & \textbf{1129} & \textbf{20.0} & \textbf{0.731} \\
 & Africa & 143 & 16.5 & 0.398 \\
 & Americas & 373 & 20.0 & 0.833 \\
 & Asia & 551 & 19.8 & 0.478 \\
 & Europe & 30 & 13.7 & 0.748 \\
 & Oceania & 32 & 36.8 & 0.672 \\
\rowcolor{lightgray} \textbf{full\_sp} & \textbf{Overall} & \textbf{1129} & \textbf{15.5} & \textbf{0.837} \\
 & Africa & 143 & 10.8 & 0.719 \\
 & Americas & 373 & 14.4 & 0.904 \\
 & Asia & 551 & 17.4 & 0.590 \\
 & Europe & 30 & 12.4 & 0.755 \\
 & Oceania & 32 & 14.8 & 0.828 \\
\rowcolor{lightgray} \textbf{full\_interact} & \textbf{Overall} & \textbf{1129} & \textbf{19.4} & \textbf{0.747} \\
 & Africa & 143 & 14.8 & 0.562 \\
 & Americas & 373 & 18.7 & 0.862 \\
 & Asia & 551 & 19.7 & 0.477 \\
 & Europe & 30 & 13.6 & 0.728 \\
 & Oceania & 32 & 36.6 & 0.649 \\
\rowcolor{lightgray} \textbf{full\_interact\_sp} & \textbf{Overall} & \textbf{1129} & \textbf{15.5} & \textbf{0.839} \\
 & Africa & 143 & 10.8 & 0.714 \\
 & Americas & 373 & 14.1 & 0.908 \\
 & Asia & 551 & 17.4 & 0.589 \\
 & Europe & 30 & 12.0 & 0.771 \\
 & Oceania & 32 & 15.0 & 0.830 \\
\rowcolor{lightgray} \textbf{elevation\_only\_sp} & \textbf{Overall} & \textbf{1129} & \textbf{15.7} & \textbf{0.834} \\
 & Africa & 143 & 11.4 & 0.688 \\
 & Americas & 373 & 14.5 & 0.903 \\
 & Asia & 551 & 17.5 & 0.584 \\
 & Europe & 30 & 12.4 & 0.757 \\
 & Oceania & 32 & 15.1 & 0.825 \\
\rowcolor{lightgray} \textbf{elevation\_c4\_sp} & \textbf{Overall} & \textbf{1129} & \textbf{15.7} & \textbf{0.834} \\
 & Africa & 143 & 11.4 & 0.685 \\
 & Americas & 373 & 14.5 & 0.903 \\
 & Asia & 551 & 17.5 & 0.585 \\
 & Europe & 30 & 12.4 & 0.758 \\
 & Oceania & 32 & 15.0 & 0.828 \\
\rowcolor{lightgray} \textbf{c4\_only\_sp} & \textbf{Overall} & \textbf{1129} & \textbf{15.9} & \textbf{0.829} \\
 & Africa & 143 & 12.2 & 0.636 \\
 & Americas & 373 & 14.7 & 0.900 \\
 & Asia & 551 & 17.7 & 0.575 \\
 & Europe & 30 & 12.4 & 0.756 \\
 & Oceania & 32 & 14.8 & 0.835 \\
\rowcolor{lightgray} \textbf{elevation\_c4\_interact\_sp} & \textbf{Overall} & \textbf{1129} & \textbf{15.7} & \textbf{0.835} \\
 & Africa & 143 & 11.5 & 0.678 \\
 & Americas & 373 & 14.4 & 0.904 \\
 & Asia & 551 & 17.5 & 0.585 \\
 & Europe & 30 & 12.4 & 0.757 \\
 & Oceania & 32 & 14.7 & 0.836 \\
"""

tab_tex.append("\n\\clearpage\n\\setcounter{table}{2}\n\\input{../tables/Table_S3_compilation_sources.tex}\n")

final = preamble + front + body_tex + "".join(tab_tex) + "".join(fig_tex) + postamble
Path("manuscript/latex/supplement.tex").write_text(final)
print(f"wrote supplement.tex ({len(final)} chars)")
print(f"figs: {sorted(fig_caps.keys())}, tabs: {sorted(tab_caps.keys())}")
print(f"unmatched cites ({len(set(UNMATCHED))}):")
for u in sorted(set(UNMATCHED))[:10]:
    print(f"  - {u}")
