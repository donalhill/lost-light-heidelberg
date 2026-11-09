"""The one results file: every number the documents state, by group, with the file it came from. Groups: design
(config/analysis.json), facts (heidelberg/facts.json, each cited), derived (heidelberg/derived.json), and one group per
results file in SOURCES. A results file's keys become the group's names, "." written "p" and other characters outside
letters, digits, and "_" written "_" (power.0.10_kappa1 becomes power.0p10_kappa1); a dict keyed by numbers also gives its
keys as a list (coverage_keys). A file marked "*" (or "**") holds several groups at its top level; later files win.
python tools/results.py OUT.json [GROUP=FILE ...]   (GROUP=FILE: another file for a group, e.g. the rebuild's calibration)"""
import json
import re
import sys

N = "heidelberg/checks"                            # this repository's numbers (Snakefile)
SOURCES = {"sample": f"{N}/S5_sample.json", "pilot": f"{N}/P1_numbers.json", "clock": f"{N}/N5_clock.json",
           "calibration": f"{N}/J3_numbers.json", "costs": f"{N}/costs.json",
           "sunlit": f"{N}/J6_sunlit.json", "forecast": f"{N}/forecast.json", "**": f"{N}/run_numbers.json"}   # "**": a file holding several groups


def leaves(x, path=""):
    if isinstance(x, dict):
        if path and x and all(re.fullmatch(r"-?\d+(\.\d+)?", str(k)) for k in x):
            yield f"{path}_keys", sorted(float(k) for k in x)
        for k, v in x.items():
            if k != "phrases":
                k = re.sub(r"\W", "_", str(k).replace(".", "p"))
                yield from leaves(v, f"{path}.{k}" if path else k)
    else:
        yield path, x


def build(over):
    src = dict(SOURCES, **over); out, frm = {}, {}
    for f, g in (("config/analysis.json", "design"), ("heidelberg/facts.json", "facts"), ("heidelberg/derived.json", "derived")):
        for sec in json.load(open(f)).items():
            if sec[0].startswith("_"):
                continue
            entries = sec[1] if g == "design" else {sec[0]: sec[1]}
            for k, e in entries.items():
                out[f"{g}.{k}"] = e["value"]; frm[f"{g}.{k}"] = f
                for x, v in e.items():                       # a number kept beside the value (an error, se)
                    if x not in ("value", "unit", "source", "reason") and isinstance(v, (int, float)):
                        out[f"{g}.{k}_{x}"] = v; frm[f"{g}.{k}_{x}"] = f
    for g, f in src.items():
        for k, v in leaves(json.load(open(f)), "" if g.startswith("*") else g):
            out[k] = v; frm[k] = f
    return out, frm


if __name__ == "__main__":
    out, frm = build(dict(a.split("=", 1) for a in sys.argv[2:]))
    json.dump(dict(values=out, sources=frm), open(sys.argv[1], "w"), indent=1, default=float)
    print(f"{len(out)} numbers in {sys.argv[1]}")
