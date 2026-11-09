"""Check 4: the shadow stage's outputs from this repository (NEW) against the previous repository's (REF), run for
run: every JSON value equal, except the cut items, which are absent from NEW (likelihood.interval).
python tools/same_shadow.py NEW_DIR REF_DIR"""
import glob
import json
import os
import sys

CUT = {("likelihood", "interval")}
new, ref = sys.argv[1:3]; bad = []; n = 0
for f in sorted(glob.glob(f"{ref}/*/*.json")):
    rel = os.path.relpath(f, ref)
    if rel.endswith("config.json"):
        continue
    g = os.path.join(new, rel)
    if not os.path.exists(g):
        bad.append(f"missing: {rel}"); continue
    a, b = json.load(open(g)), json.load(open(f))
    for top, key in CUT:
        if key in a.get(top, {}):
            bad.append(f"{rel}: {top}.{key} should be absent")
        b.get(top, {}).pop(key, None)
    n += 1
    if a != b:
        diff = [k for k in set(a) | set(b) if a.get(k) != b.get(k)]
        bad.append(f"{rel}: differs in {sorted(diff)}")
for x in bad:
    print(x)
print(f"{n} outputs compared: " + ("equal but for the cut items" if not bad else f"{len(bad)} differences"))
sys.exit(1 if bad else 0)
