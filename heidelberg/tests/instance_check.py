"""The primary selection made here (offline from the dry run's frozen answers) against the dry run's results on the
instance, whose code also decided on AllWISE (removed since, intended difference 14): the sealed candidates row by row
(every current cut's flag, the remote evaluation, positions; a row passes here exactly when it passed there or failed
there on AllWISE alone), the sealed candidate counts, and the open per-pair records (the injections; per cut, the chance
rows failing it; the pairing cut's row, which must be empty for chance rows).
python -E -s -m heidelberg.tests.instance_check NEW_SELECT_DIR DRY_SEALED_DIR DRY_OPEN_SELECT_DIR"""
import glob, json, os, sys
import numpy as np, pandas as pd
NEW, DRY_SEALED, DRY_OPEN = sys.argv[1:4]
from heidelberg.lib.selection import CUTS
bad = []; nrow = 0; fin_old = fin_new = 0
for f in sorted(glob.glob(f"{DRY_OPEN}/*_select.json")):
    o = json.load(open(f)); p = o["pair"]; n = json.load(open(f"{NEW}/{p}_select.json"))
    if "refused" in o or "refused" in n:
        if ("refused" in o) != ("refused" in n): bad.append(f"{p}: refusal differs")
        continue
    if o.get("injected") != n.get("injected"): bad.append(f"{p}: injected {o.get('injected')} vs {n.get('injected')}")   # the chance estimate gains the rows AllWISE alone removed
    for pop in ("chance",):
        a, b = o["cutflow"].get(pop, {}), n["cutflow"].get(pop, {})
        if b.get("pairing", [0, 0]) != [0, 0]: bad.append(f"{p}: chance rows fail the pairing cut")   # it judges injections only
        b = {k: v for k, v in b.items() if k != "pairing"}                                   # the row the pairing cut adds (difference 8)
        if any(a.get(c, [0])[0] != b.get(c, [0])[0] for c in CUTS if c != "pairing") or a.get("rows") != b.get("rows"):
            bad.append(f"{p}: cut-flow {pop} differs")
        if b.get("passed") != a.get("passed") + a.get("allwise", [0, 0])[1]: bad.append(f"{p}: chance passes {a.get('passed')} + AllWISE alone vs {b.get('passed')}")
    so = json.load(open(f"{DRY_SEALED}/{p}_select.json")) if os.path.exists(f"{DRY_SEALED}/{p}_select.json") else {}
    sn = json.load(open(f"{NEW}/sealed/{p}_select.json"))
    if so.get("candidates") != sn.get("candidates"): bad.append(f"{p}: sealed candidates {so.get('candidates')} vs {sn.get('candidates')}")
    fin_old += so.get("final", 0); fin_new += sn.get("final", 0)
    A = pd.read_parquet(f"{DRY_SEALED}/{p}_candidates.parquet"); B = pd.read_parquet(f"{NEW}/sealed/{p}_candidates.parquet")
    if len(A) == len(B) and not np.array_equal(B["passed"].to_numpy(bool) & ~A["fail_allwise"].to_numpy(bool), A["passed"].to_numpy(bool)):
        bad.append(f"{p}: a verdict differs beyond AllWISE")
    if sn.get("final") != int(B["passed"].sum()): bad.append(f"{p}: sealed final {sn.get('final')} vs {int(B['passed'].sum())} rows")
    cols = [f"fail_{c}" for c in CUTS if f"fail_{c}" in A] + ["remote_evaluated"]
    if len(A) != len(B) or not A[cols].reset_index(drop=True).equals(B[cols].reset_index(drop=True)):
        d = {c: int((A[c].to_numpy() != B[c].to_numpy()).sum()) for c in cols if len(A) == len(B) and (A[c].to_numpy() != B[c].to_numpy()).any()}
        bad.append(f"{p}: candidate flags differ {d if d else (len(A), len(B))}")
    for c in ("ra", "dec", "ra_bo", "dec_bo"):
        if len(A) == len(B) and not np.array_equal(A[c].to_numpy(), B[c].to_numpy()): bad.append(f"{p}: {c} differs")
    nrow += len(A)
print(f"{nrow} candidate rows compared; finals {fin_old} (instance) vs {fin_new} (laptop)")
for b in bad[:40]: print("DIFFERS", b)
print("IDENTICAL to the instance" if not bad else f"{len(bad)} differences"); sys.exit(1 if bad else 0)
