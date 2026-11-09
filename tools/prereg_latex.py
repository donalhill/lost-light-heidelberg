"""The rendered pre-registration (draft/preregistration_v2.md) as a LaTeX document, for a tidy PDF: headings, paragraphs,
lists, tables (booktabs), the model's formula displayed, in Computer Modern (Latin Modern), and the symbols set as mathematics. Compile with tectonic.
python3 tools/prereg_latex.py draft/preregistration_v2.md OUT.tex [DATE_LINE]"""
import re
import sys

SRC, OUT = sys.argv[1:3]
DATE = sys.argv[3] if len(sys.argv) > 3 else ""

MATH = [  # (text as rendered, LaTeX), longest first; applied after escaping
    ("λᵢ(x) = Bᵢ wᵢ(x) + η Eᵢ [1 − κ s(x)].", None),               # the displayed formula, handled on its own
    ("η / (η + B)", r"$\eta/(\eta + B)$"), ("1 − κs", r"$1 - \kappa s$"), ("κ > 0", r"$\kappa > 0$"),
    ("κ = 0.5", r"$\kappa = 0.5$"), ("κ = 1", r"$\kappa = 1$"), ("κ = 0", r"$\kappa = 0$"),
    ("s = 1/3", r"$s = 1/3$"), ("s = 1", r"$s = 1$"), ("s = 0", r"$s = 0$"),
    ("H0:", r"$H_0$:"), ("H1 (primary):", r"$H_1$ (primary):"), ("shadow share s,", r"shadow share $s$,"),
    ("shadow share, s,", r"shadow share, $s$,"), ("We compute s where", r"We compute $s$ where"), ("and s then counts", r"and $s$ then counts"),
    ("In this expression, B is", r"In this expression, $B$ is"), ("B converts", r"$B$ converts"), ("use B =", r"use $B =$"), ("statistic, q,", r"statistic, $q$,"), ("statistic, D,", r"statistic, $D$,"),
    ("pair's n selected", r"pair's $n$ selected"), ("average to n times", r"average to $n$ times"),
    ("position x on pair i", r"position $x$ on pair $i$"), ("Each plate pair i has", r"Each plate pair $i$ has"),
    ("(config/analysis.json)", r"(\texttt{config/analysis.json})"),
    ("H0", r"$H_0$"), ("H1", r"$H_1$"), ("H2", r"$H_2$"), ("Bᵢ", r"$B_i$"), ("wᵢ", r"$w_i$"), ("Eᵢ", r"$E_i$"), ("κ", r"$\kappa$"), ("η", r"$\eta$"),
]


def esc(t):
    t = t.replace("\\", r"\textbackslash{}")
    for a, b in (("&", r"\&"), ("%", r"\%"), ("$", r"\$"), ("#", r"\#"), ("_", r"\_"), ("{", r"\{"), ("}", r"\}"),
                 ("~", r"\textasciitilde{}"), ("^", r"\textasciicircum{}")):
        t = t.replace(a, b)
    return t


def inline(t):
    t = esc(t); keep = []
    for a, b in MATH:                                                       # placeholders first, so no rule rematches
        if b and esc(a) in t:
            keep.append(b); t = t.replace(esc(a), f"\x00{len(keep) - 1}\x00")
    t = re.sub(r'"([^"]+)"', r"``\1''", t)                                  # quotation marks
    t = re.sub(r"(\d) x (\d)", r"\1\\,$\\times$\\,\2", t)                   # 6 x 6 degrees
    t = t.replace("−", r"\textminus{}")
    t = re.sub(r"\b(et al|e\.g|i\.e|cf|Dr)\. ", lambda m: m.group(1) + ".\\ ", t)   # an abbreviation's stop is not a sentence's
    t = re.sub(r"\*\*(.+?)\*\*", r"\\textbf{\1}", t)
    t = re.sub(r"(?<![\w*])\*(?!\s)(.+?)(?<!\s)\*(?![\w*])", r"\\emph{\1}", t)        # *journal*
    return re.sub("\x00(\\d+)\x00", lambda m: keep[int(m.group(1))], t)


def table(lines):
    rows = [[c.strip() for c in ln.strip().strip("|").split("|")] for ln in lines if not re.match(r"^\|\s*-", ln)]
    n = len(rows[0])
    L, Rt = r">{\raggedright\arraybackslash}X", r">{\raggedleft\arraybackslash}X"
    spec = (L + Rt * (n - 1)) if n > 2 else (L + r">{\raggedleft\arraybackslash}p{2cm}")
    out = [r"\begin{center}\small", rf"\begin{{tabularx}}{{\linewidth}}{{{spec}}}", r"\toprule",
           " & ".join(rf"\textbf{{{inline(c)}}}" for c in rows[0]) + r" \\", r"\midrule"]
    for r_ in rows[1:]:
        last = r_[0].lower().startswith("total")
        out.append(" & ".join((rf"\textbf{{{inline(c)}}}" if last else inline(c)) for c in r_) + r" \\")
        if last:
            out.insert(-1, r"\midrule")
    return out + [r"\bottomrule", r"\end{tabularx}", r"\end{center}"]


text = re.sub(r"<!--.*?-->", "", open(SRC).read(), flags=re.S).strip()
blocks = re.split(r"\n\s*\n", text)
title = blocks[0].lstrip("# ").strip()
body = []
in_refs = False
for b in blocks[1:]:
    lines = [ln for ln in b.split("\n") if ln.strip()]
    if not lines:
        continue
    first = lines[0]
    if first.startswith("## "):
        h = first[3:].strip(); in_refs = h.lower().endswith("references")
        body.append(rf"\section*{{{inline(h)}}}")
        lines = lines[1:]
        if not lines:
            continue
        first = lines[0]
    if first.startswith("#### "):
        body.append(rf"\subsubsection*{{{inline(first[5:].strip())}}}"); lines = lines[1:]
        if not lines:
            continue
        first = lines[0]
    if first.startswith("### "):
        body.append(rf"\subsection*{{{inline(first[4:].strip())}}}"); lines = lines[1:]
        if not lines:
            continue
        first = lines[0]
    if first.startswith("|"):
        body += table(lines)
    elif first.startswith("λᵢ(x)"):
        body.append(r"\[ \lambda_i(x) = B_i\, w_i(x) + \eta\, E_i\, \bigl[1 - \kappa\, s(x)\bigr]. \]")
    elif re.match(r"^\d+\. ", first):                                       # the numbered references
        body.append(r"\begin{list}{}{\leftmargin=2em\labelwidth=1.6em\labelsep=0.4em\itemsep=2pt\parsep=0pt}")
        for ln in lines:
            m = re.match(r"^(\d+)\. (.*)$", ln)
            body.append(rf"\item[{{[{m.group(1)}]}}] {inline(m.group(2))}")
        body.append(r"\end{list}")
    elif first.startswith("- "):
        env = r"\begin{list}{}{\leftmargin=1.2em\itemindent=-1.2em\itemsep=2pt\parsep=0pt}" if in_refs else r"\begin{itemize}\setlength\itemsep{2pt}"
        body.append(env)
        items = []
        for ln in lines:
            if ln.startswith("- "):
                items.append(ln[2:])
            else:
                items[-1] += " " + ln.strip()
        body += [rf"\item {inline(i)}" for i in items]
        body.append(r"\end{list}" if in_refs else r"\end{itemize}")
    else:
        p = " ".join(ln.strip() for ln in lines)
        p = re.sub(r"^\((\d)\) ([^.]+)\.", lambda m: rf"\textbf{{({m.group(1)}) {m.group(2)}.}}", p)   # the abstract's prompts
        body.append(inline(p) if not p.startswith(r"\textbf") else p[:p.index("}") + 1] + inline(p[p.index("}") + 1:]))
    body.append("")

doc = r"""\documentclass[11pt,a4paper]{article}
\usepackage[margin=24mm]{geometry}
\usepackage{fontspec}
\usepackage{unicode-math}
\setmainfont{lmroman10}[Extension=.otf, UprightFont=*-regular, BoldFont=*-bold, ItalicFont=*-italic, BoldItalicFont=*-bolditalic]
\setmathfont{latinmodern-math.otf}
\usepackage{microtype}
\usepackage{booktabs}
\usepackage{tabularx}
\usepackage{array}
\emergencystretch=2em
\usepackage[hidelinks]{hyperref}
\usepackage{titlesec}
\titleformat{\section}{\large\bfseries}{}{0pt}{}
\titleformat{\subsection}{\normalsize\bfseries}{}{0pt}{}
\titleformat{\subsubsection}{\normalsize\itshape}{}{0pt}{}
\titlespacing*{\section}{0pt}{14pt}{6pt}
\titlespacing*{\subsection}{0pt}{10pt}{4pt}
\setlength{\parindent}{0pt}
\setlength{\parskip}{5pt}
\linespread{1.08}
\begin{document}
\begin{center}
{\LARGE\bfseries """ + inline(title) + r"""\par}
\vspace{6pt}
{\large Donal Hill\par}
""" + (rf"\vspace{{2pt}}{{\small {esc(DATE)}\par}}" if DATE else "") + r"""
\end{center}
\vspace{4pt}
""" + "\n".join(body) + r"""
\end{document}
"""
open(OUT, "w").write(doc)
print(f"{OUT}: {len(body)} lines")
