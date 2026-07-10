"""Convert leafwax_gca_text_2026-05-05.md → elsarticle main.tex (v2 with all Unicode/escape fixes)."""
import re, json, sys, unicodedata
from pathlib import Path

KEYS = set(json.loads(Path("manuscript/latex/_keys.json").read_text()))
UNMATCHED = []

def normalize_surname(s):
    s = unicodedata.normalize("NFKD", s).encode("ascii", "ignore").decode("ascii").lower()
    return re.sub(r"[^a-z]", "", s)

def lookup_key(surname, year_letter):
    n = normalize_surname(surname)
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

def md_inline(s):
    # Markdown-escaped chars FIRST, before they get interpreted
    s = s.replace(r"\'", "'")
    s = s.replace(r"\<", "<").replace(r"\>", ">").replace(r"\|", "|")
    s = s.replace(r"\~", "~")  # markdown literal tilde — but we'll convert ~X~ subscripts later
    s = s.replace(r"\*", "*")  # markdown literal star (rare in this doc)
    # ε~app~ specific common form
    s = s.replace("ε~app~", r"$\varepsilon_{\mathrm{app}}$")
    # Combining circumflex BEFORE we tokenize: R̂ → $\hat{R}$ etc.
    # `̂` is U+0302 combining circumflex; preceding char is the base
    s = re.sub(r"([A-Za-z])̂", r"$\\hat{\1}$", s)
    # superscript markdown ^X^
    s = re.sub(r"\^([^\s^]+?)\^", r"\\textsuperscript{\1}", s)
    # subscript markdown ~X~
    s = re.sub(r"~([^\s~]+?)~", r"\\textsubscript{\1}", s)
    # bold + italic (after sub/sup so the bold/emph wraps text not markup)
    s = re.sub(r"\*\*(.+?)\*\*", r"\\textbf{\1}", s)
    s = re.sub(r"(?<!\\)\*([^*\n]+?)\*", r"\\emph{\1}", s)
    repl = {
        # Greek
        "δ": r"$\delta$", "Δ": r"$\Delta$", "ε": r"$\varepsilon$",
        "α": r"$\alpha$", "β": r"$\beta$", "λ": r"$\lambda$",
        "ρ": r"$\rho$", "σ": r"$\sigma$", "μ": r"$\mu$", "τ": r"$\tau$",
        "χ": r"$\chi$", "ν": r"$\nu$", "η": r"$\eta$", "θ": r"$\theta$",
        # ‰ deg etc.
        "‰": r"\textperthousand{}",
        # subscript/superscript Unicode digits
        "⁰": r"$^{0}$", "¹": r"$^{1}$", "²": r"$^{2}$", "³": r"$^{3}$",
        "⁴": r"$^{4}$", "⁵": r"$^{5}$", "⁶": r"$^{6}$", "⁷": r"$^{7}$",
        "⁸": r"$^{8}$", "⁹": r"$^{9}$",
        "₀": r"$_{0}$", "₁": r"$_{1}$", "₂": r"$_{2}$", "₃": r"$_{3}$",
        "₄": r"$_{4}$", "₅": r"$_{5}$", "₆": r"$_{6}$", "₇": r"$_{7}$",
        "₈": r"$_{8}$", "₉": r"$_{9}$",
        "ᵢ": r"$_{i}$", "ⱼ": r"$_{j}$",
        "ₖ": r"$_{k}$", "ₙ": r"$_{n}$",
        # punctuation
        "−": "-", "–": "--", "—": "---",
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
    }
    for k, v in repl.items():
        s = s.replace(k, v)
    # Tilde-as-approximation in source: "\~22%" (md-escaped) → already "~22%" after \~ replace
    # Now interpret remaining ~ as approximation if followed by digit/word
    s = re.sub(r"(?<!\$)~(?=[\d])", r"$\\sim$", s)
    s = re.sub(r"(?<!\$)~(?=[A-Za-z])", r"$\\sim$", s)  # only in non-subscript context — sub already converted
    # Plain `<` and `>` outside math: convert to \textless / \textgreater (BUT not inside <abradley@..>)
    # First handle email: <abradley@wustl.edu> -> \href{mailto:...}{...}
    s = re.sub(r"<([A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,})>", r"\\href{mailto:\1}{\1}", s)
    # Bare `<` and `>` to math mode
    s = re.sub(r"(?<![\\$])<", r"$<$", s)
    s = re.sub(r"(?<![\\$])>", r"$>$", s)
    # Pipe in absolute value: |x|
    # Wrap |...|>0.7 → $|...|>0.7$ — too risky to autodetect; leave |s as-is but escape if outside math.
    # Final: escape underscores OUTSIDE math mode. Walk through string char-by-char tracking $ pairs.
    out = []
    in_math = False
    i = 0
    while i < len(s):
        c = s[i]
        if c == '$':
            in_math = not in_math
            out.append(c)
            i += 1
            continue
        if c == '_' and not in_math:
            # Don't escape if already escaped (\_)
            if i > 0 and s[i-1] == '\\':
                out.append(c); i += 1; continue
            out.append(r'\_'); i += 1; continue
        if c == '#' and not in_math:
            # Escape standalone # (but not in URLs etc.) — be conservative
            out.append(r'\#'); i += 1; continue
        if c == '&' and not in_math:
            # Escape & — but it might already be escaped \&
            if i > 0 and s[i-1] == '\\':
                out.append(c); i += 1; continue
            out.append(r'\&'); i += 1; continue
        out.append(c); i += 1
    return "".join(out)

def convert_text(text):
    text = replace_narrative_cites(text)
    text = replace_paren_cites(text)
    text = md_inline(text)
    return text

SRC = Path("manuscript/drafts/leafwax_gca_text_2026-05-05.md").read_text()
body_text = SRC.split("**References**", 1)[0]

def scan_body(body):
    lines = body.splitlines()
    out, i, in_eq, eq_buf = [], 0, False, []
    HEADER = re.compile(r"^\*\*(\d+(?:\.\d+)*)\.?\s+([^*]+?)\*\*\s*$")
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
        eql_m = re.match(r"^\(\s*Equation\s+(\d+)\s*\)\s*$", line.strip())
        if eql_m:
            out.append(("equation_label", eql_m.group(1))); i += 1; continue
        sp = SPECIAL.match(line.strip())
        if sp:
            out.append(("special", sp.group(1).strip().rstrip(":")))
            i += 1; continue
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

abs_idx = key_idx = first_header_idx = None
for k, el in enumerate(grouped):
    if el[0] == "special":
        nm = el[1].lower()
        if nm == "abstract" and abs_idx is None: abs_idx = k
        elif nm in ("keywords", "keyword") and key_idx is None: key_idx = k
    if el[0] == "header" and first_header_idx is None: first_header_idx = k

abstract_text = ""
if abs_idx is not None:
    paras = []
    for el in grouped[abs_idx+1:]:
        if el[0] in ("special", "header"): break
        if el[0] == "para":
            # Strip embedded "**Keywords:**..." line  
            t = el[1]
            # If a paragraph starts with "**Keywords:**" treat as keyword line, skip
            if re.match(r"^\*\*Keywords?:\*\*", t.strip()):
                continue
            paras.append(t)
    abstract_text = "\n\n".join(paras)

keywords_text = "leaf wax, hydrogen isotopes, paleoclimate, spatial autocorrelation"
# Override from .md if a Keywords line is found
m = re.search(r"\*\*Keywords?:\*\*\s*([^\n]+)", SRC)
if m: keywords_text = m.group(1).strip()

title_text = ""
m = re.match(r"^\*\*(.+?)\*\*\s*$", SRC.splitlines()[0])
if m: title_text = m.group(1)

body_chunks = []
sections_levels = {1: r"\section", 2: r"\subsection", 3: r"\subsubsection"}
for el in grouped[first_header_idx:]:
    t = el[0]
    if t == "header":
        depth, title, num = el[1], el[2], el[3]
        cmd = sections_levels.get(depth, r"\subsubsection")
        body_chunks.append(f"\n{cmd}{{{md_inline(title)}}}\n\\label{{sec:s{num.replace('.','_')}}}\n")
    elif t == "special":
        nm = el[1].lower()
        if nm in ("abstract", "keywords", "keyword"): continue
        if nm == "acknowledgments": body_chunks.append("\n\\section*{Acknowledgments}\n")
        elif nm == "data availability": body_chunks.append("\n\\section*{Data Availability}\n")
        elif nm == "appendix a: supplementary materials": body_chunks.append("\n\\section*{Appendix A. Supplementary Materials}\n")
        elif nm in ("figure captions", "table titles"): break
        else: body_chunks.append(f"\n\\paragraph{{{el[1]}}}\n")
    elif t == "para":
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
abstract_tex = convert_text(abstract_text or "")
keywords_tex_raw = convert_text(keywords_text or "")
# Convert commas to \sep for elsarticle keyword block
keywords_tex = " \\sep ".join([k.strip() for k in keywords_tex_raw.split(",") if k.strip()])
title_tex = convert_text(title_text or "")

# --- Figure captions (handle the ** in titles like Fig 3) ---
fig_caps = {}
fig_block_match = re.search(r"\*\*Figure Captions\*\*\s*\n(.*?)(?=\n\*\*Table titles|\n\*\*References|\Z)", SRC, re.DOTALL)
if fig_block_match:
    fc = fig_block_match.group(1)
    # Split by paragraph boundary preceding "**Figure N."
    blocks = re.split(r"\n\n(?=\*\*Figure \d+\.)", fc.strip())
    for b in blocks:
        m = re.match(r"\*\*Figure (\d+)\.\s+(.*)", b, re.DOTALL)
        if not m: continue
        n = int(m.group(1))
        rest = m.group(2)
        # Find closing ** that ends the bold "short title" — it's the FIRST `**` that is NOT part of a nested **bold**
        # Heuristic: scan for `**` and count depth
        # Simpler: short title typically ends at `**.` (period after closing) at top-level. Search for `\*\*\s` followed by sentence text
        # Robust approach: find the LAST `**` followed by a space/period+space at the top level — use depth tracker
        depth = 0
        end_short = None
        i = 0
        while i < len(rest):
            if rest[i:i+2] == "**":
                if depth == 0:
                    depth = 1
                else:
                    depth = 0
                    if i+2 < len(rest) and rest[i+2] in ".\n " and depth == 0:
                        end_short = i+2
                        break
                i += 2
                continue
            i += 1
        if end_short is None:
            # Fallback: split on first ".**\s"
            mm = re.search(r"\.\*\*\s", rest)
            if mm:
                end_short = mm.end() - 1
            else:
                end_short = 80
        short = rest[:end_short].strip().rstrip("*").rstrip(".").strip()
        full = rest[end_short:].lstrip(". \n*").strip()
        fig_caps[n] = (short, full)

tab_titles = {}
tt_match = re.search(r"\*\*Table titles\*\*\s*\n(.*?)(?=\n\*\*References|\Z)", SRC, re.DOTALL)
if tt_match:
    for m in re.finditer(r"\*\*Table (\d+)\.\*\*\s+([^\n]+)", tt_match.group(1)):
        tab_titles[int(m.group(1))] = m.group(2).strip().rstrip(".")

preamble = r"""\documentclass[5p,times,authoryear]{elsarticle}

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

% GCA wants continuous line numbers
\modulolinenumbers[1]
\linenumbers

\graphicspath{{../figures/main_figs/}}

\journal{Geochimica et Cosmochimica Acta}
"""

front = f"""
\\begin{{document}}

\\begin{{frontmatter}}

\\title{{{title_tex}}}

\\author[wustl]{{Alexander S. Bradley\\corref{{cor1}}}}
\\ead{{abradley@wustl.edu}}
\\cortext[cor1]{{Corresponding author.}}

\\address[wustl]{{Department of Earth, Environmental, and Planetary Sciences,
                Washington University in St. Louis,
                1 Brookings Drive, Saint Louis, Missouri 63130, USA}}

\\begin{{abstract}}
{abstract_tex}
\\end{{abstract}}

\\begin{{keyword}}
{keywords_tex}
\\end{{keyword}}

\\end{{frontmatter}}
"""

postamble = r"""

\bibliographystyle{elsarticle-harv}
\bibliography{references}

\end{document}
"""

fig_files = {
    1: "Figure_01_combined.pdf",
    2: "Figure_02_all_environmental_variables.png",
    3: "figure_03_spatial_confounding.pdf",
    4: "Figure_04.pdf",
    5: "Figure_05_detection_thresholds.png",
}
fig_tex = ["\n\\clearpage\n"]
for n in sorted(fig_caps.keys()):
    short, full = fig_caps[n]
    fname = fig_files.get(n, f"Figure_{n:02d}.pdf")
    short_clean = convert_text(short)
    full_clean = convert_text(full)
    fig_tex.append(rf"""
\begin{{figure}}[p]
  \centering
  \includegraphics[width=0.95\linewidth]{{{fname}}}
  \caption[{short_clean}]{{\textbf{{{short_clean}.}} {full_clean}}}
  \label{{fig:fig{n}}}
\end{{figure}}
""")

tab_tex = ["\n\\clearpage\n"]
tab_inputs = {
    1: "../tables/table1_model_performance.tex",
    2: None,
    3: "../tables/table3_variance_decomp.tex",
    4: "../tables/table4_environmental.tex",
}
for n in sorted(tab_titles.keys()):
    title = tab_titles[n]
    title_clean = convert_text(title)
    src = tab_inputs.get(n)
    if src:
        tab_tex.append(f"\n\\input{{{src}}}\n")
    else:
        # Table 2 uses a complex longtable structure pulled from the standalone wrapper
        tab_tex.append(rf"""
\clearpage
\begin{{landscape}}
\small
\begin{{longtable}}{{@{{}}l>{{\raggedleft\arraybackslash}}p{{2.2cm}}>{{\raggedleft\arraybackslash}}p{{2.2cm}}>{{\raggedleft\arraybackslash}}p{{2.3cm}}>{{\raggedleft\arraybackslash}}p{{2.2cm}}>{{\raggedleft\arraybackslash}}p{{1.8cm}}>{{\raggedleft\arraybackslash}}p{{2.2cm}}cc@{{}}}}
\caption{{\textbf{{Table {n}.}} {title_clean}}}\label{{tab:tab{n}}} \\
\toprule
\textbf{{Model}} & \textbf{{RMSE}} & \textbf{{R}}$^2$ & \textbf{{$\beta_0$}} & \textbf{{$\beta_{{\text{{OIPC}}}}$}} & \textbf{{$\lambda_{{\text{{int}}}}$}} & \textbf{{GP scale}} & \textbf{{Knot}} & \textbf{{Knot}} \\
 & \textbf{{(\textperthousand)}} & & \textbf{{(\textperthousand)}} & & \textbf{{(km)}} & \textbf{{(km)}} & \textbf{{Slope SD}} & \textbf{{Int. SD}} \\
 & & & & & & & & \textbf{{(\textperthousand)}} \\
\midrule
\endfirsthead
\multicolumn{{9}}{{c}}{{Table {n} (continued)}} \\
\toprule
\textbf{{Model}} & \textbf{{RMSE}} & \textbf{{R}}$^2$ & \textbf{{$\beta_0$}} & \textbf{{$\beta_{{\text{{OIPC}}}}$}} & \textbf{{$\lambda_{{\text{{int}}}}$}} & \textbf{{GP scale}} & \textbf{{Knot}} & \textbf{{Knot}} \\
 & \textbf{{(\textperthousand)}} & & \textbf{{(\textperthousand)}} & & \textbf{{(km)}} & \textbf{{(km)}} & \textbf{{Slope SD}} & \textbf{{Int. SD}} \\
 & & & & & & & & \textbf{{(\textperthousand)}} \\
\midrule
\endhead
\bottomrule
\endlastfoot
\input{{../tables/table2_global_params_body.tex}}
\end{{longtable}}
\end{{landscape}}
""")

final = preamble + front + body_tex + "".join(tab_tex) + "".join(fig_tex) + postamble
Path("manuscript/latex/main.tex").write_text(final)
print(f"wrote main.tex ({len(final)} chars)")
print(f"figs found: {sorted(fig_caps.keys())}, tabs: {sorted(tab_titles.keys())}")
print(f"unmatched cites: {len(set(UNMATCHED))}")
