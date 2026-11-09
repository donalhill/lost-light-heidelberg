"""A document rendered from its template and the results file (tools/results.py). A placeholder is {{key}} or
{{key|format}}: key a number's name, or an element of a list as key.N (key.N.M in a list of lists); formats, chained with
commas: int (2,768), words (up to twelve: six), dN (N decimals), pctN (a percentage with N decimals), neg (a typographic
minus), cap (the first letter in capitals, after words: Five), list ("a, b, and c"), len, none ("none", and the render stops unless the value is zero). Rounding is half up from the shortest decimal form (0.145 to 0.15). An unknown key or
format stops the render.
python tools/render.py TEMPLATE.md RESULTS.json OUT.md"""
import json
import re
import sys
from decimal import ROUND_HALF_UP, Decimal

WORDS = "zero one two three four five six seven eight nine ten eleven twelve".split()
PH = re.compile(r"\{\{\s*([\w.]+)\s*(?:\|\s*([\w,]+)\s*)?\}\}")


def half_up(v, n=0):
    return (v if isinstance(v, Decimal) else Decimal(repr(float(v)))).quantize(Decimal(1).scaleb(-n), rounding=ROUND_HALF_UP)


def lookup(R, key):
    if key in R:
        return R[key]
    parts = key.split(".")
    for n in range(len(parts) - 1, 0, -1):
        base, idx = ".".join(parts[:n]), parts[n:]
        if base in R and all(p.isdigit() for p in idx):
            v = R[base]
            for p in idx:
                v = v[int(p)]
            return v
    raise KeyError(f"no number named {key!r}")


def fmt(v, spec):
    s = None
    for f in spec.split(",") if spec else []:
        if f == "int":
            s = f"{int(half_up(v)):,}"
        elif f == "words":
            n = int(half_up(v)); s = WORDS[n] if 0 <= n <= 12 else f"{n:,}"
        elif re.fullmatch(r"d\d", f):
            s = f"{half_up(v, int(f[1])):.{int(f[1])}f}"
        elif re.fullmatch(r"pct\d", f):
            s = f"{half_up(Decimal(repr(float(v))) * 100, int(f[3])):.{int(f[3])}f}"
        elif f == "len":
            s = f"{len(v):,}"
        elif f == "list":
            xs = [f"{x:g}" for x in v]; s = ", ".join(xs[:-1]) + ", and " + xs[-1] if len(xs) > 2 else " and ".join(xs)
        elif f == "none":                       # a count stated as "none": the render stops if the value is not zero
            if v != 0:
                raise ValueError(f"format none: the value is {v}, not zero; the sentence needs rewording")
            s = "none"
        elif f == "cap":
            s = s if s is not None else str(v); s = s[:1].upper() + s[1:]
        elif f == "neg":
            s = (s if s is not None else str(v)).replace("-", "−")
        else:
            raise ValueError(f"unknown format {f!r}")
    return s if s is not None else (f"{v:,}" if isinstance(v, int) and not isinstance(v, bool) else str(v))


if __name__ == "__main__":
    tpl, res, out = sys.argv[1:4]
    R = json.load(open(res))["values"]; text = open(tpl).read()
    missing = []
    for m in PH.finditer(text):
        try:
            lookup(R, m.group(1))
        except (KeyError, IndexError, TypeError):
            missing.append(m.group(1))
    if missing:
        sys.exit(f"{len(missing)} placeholders with no number: " + ", ".join(sorted(set(missing))))
    rendered = PH.sub(lambda m: fmt(lookup(R, m.group(1)), m.group(2)), text)
    open(out, "w").write(rendered); print(f"{out}: {len(PH.findall(text))} placeholders filled")
