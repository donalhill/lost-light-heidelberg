"""Stage 3, per pair (laptop): the three populations of measure_pair.py's products (candidates, chance pairs, injected
glints) as pair tables, each through select() (heidelberg/lib/selection.py), with the catalogue answers from the cache
(frozen for the record run: --offline). A pair is refused, with its reason, when it is not accepted (solve), when its
injection fidelity check failed, or when its products predate the forced and shifted measurements. Writes, per pair, to
OUT_DIR: PAIR_flags.parquet (every row of every population with its fail flags), PAIR_eff.csv, and PAIR_select.json
(cut-flows, counts, and the cache files used). Cut-flow entries are [fails the cut, fails it alone]
(lib/decision.cutflow);
for the remote catalogue cuts (selection.REMOTE) the first is among the rows remote_evaluated (those failing at most one
other cut).
python -m heidelberg.stages.select PRODUCTS_DIR OUT_DIR [PAIRS...] --xcache DIR --asteroids DIR --scanner-marks CSV
    --snapshots DIR [--offline] [--blind] [--jobs N] [--ps1 hats|mast]
--blind: the record run (candidates' rows and counts written only to OUT_DIR/sealed/, which the job time-locks; nothing
per pair about candidates printed); --offline: every catalogue answer from the cache (frozen); --ps1: Pan-STARRS from
the HATS copy on AWS (default) or MAST's API. Every pair accepted by its solve must have its fidelity report
(PAIR_fidelity.json, heidelberg/stages/fidelity.py). The parameters come from config/analysis.json."""
import argparse
import glob
import json
import os

import numpy as np
import pandas as pd

from heidelberg.config import A
from heidelberg.lib import catalogues as CAT
from heidelberg.lib import decision as DEC
from heidelberg.lib import selection as SEL

ap = argparse.ArgumentParser()                     # parsed at import: the spawned workers re-import with the same arguments
ap.add_argument("products"); ap.add_argument("out"); ap.add_argument("pairs", nargs="*")
ap.add_argument("--xcache", required=True); ap.add_argument("--asteroids", required=True)
ap.add_argument("--scanner-marks", required=True); ap.add_argument("--snapshots", required=True)
ap.add_argument("--offline", action="store_true"); ap.add_argument("--blind", action="store_true")
ap.add_argument("--jobs", type=int, default=1); ap.add_argument("--ps1", choices=("hats", "mast"), default="hats")
ARGS = ap.parse_args(); PD, OD, ONLY = ARGS.products, ARGS.out, set(ARGS.pairs)
os.makedirs(OD, exist_ok=True)
SEL.configure(source="primary", xcache=ARGS.xcache, asteroids=ARGS.asteroids, scanner_marks=ARGS.scanner_marks,
              offline=ARGS.offline, rule="symmetric", snapshots=ARGS.snapshots)
CAT.PS1_SOURCE = ARGS.ps1
BLIND = ARGS.blind                                 # the record run: candidates' per-pair results only in OUT_DIR/sealed/ (time-locked)
SEALED = f"{OD}/sealed"
if BLIND:
    os.makedirs(SEALED, exist_ok=True)


def one(f):
    base = f[:-5]; pair = os.path.basename(base); o = json.load(open(f)); ob = f"{OD}/{pair}"
    za, zb = np.load(base + "_a.npz"), np.load(base + "_b.npz"); g = np.load(base + "_gaia.npz")
    refuse, fid_src = None, "PAIR_fidelity.json"     # the fidelity report of heidelberg/stages/fidelity.py
    if min(x["gaia_matched"] for x in o["solve"]) < A.solve.solve_min_stars or not all(x.get("accepted", False) for x in o["solve"]):
        refuse = "solve not accepted"
    elif not os.path.exists(base + "_fidelity.json"):
        raise RuntimeError(f"{pair}: no fidelity report (run heidelberg.stages.fidelity on the products first)")
    elif "zshift" not in za.files or "zshift" not in zb.files:
        refuse = "products predate the shifted forced measurements"
    elif not (fid := json.load(open(base + "_fidelity.json"))).get("passed"):
        refuse = f"injection fidelity: {fid}"
    if refuse:
        r = dict(pair=pair, refused=refuse); json.dump(r, open(ob + "_select.json", "w"), indent=1, default=float); return r
    SEL.USED.clear(); C = SEL.context(base, o, za, zb, g)
    tables = []
    ia, jb = SEL.candidates(za, zb); Tc = SEL.det_table(za, zb, ia, jb).assign(population="candidate", ia=ia, jb=jb); tables.append(Tc)
    for k, T_ in enumerate(SEL.chance_tables(za, zb)):
        tables.append(T_.assign(population="chance", shift=k))
    for r_ in (list(za["xshift_radii"]) if "xshift_radii" in za.files else []):     # the same estimate at other lags (a test, not the estimate)
        for k, T_ in enumerate(SEL.chance_tables(za, zb, float(r_)) or []):
            tables.append(T_.assign(population=f"chance_{int(r_)}", shift=k))
    J = pd.read_parquet(base + "_inj.parquet")
    if len(J) and "found_a" in J:
        Jf, Ti = SEL.injection_table(J); tables.append(Ti.assign(population="injection", inj=Ti.index, G=Jf.G, dur=Jf.dur, sharpen=Jf.sharpen))
    T = pd.concat(tables, ignore_index=True); F = SEL.select(T, C); ev = F.pop("remote_evaluated")
    T = pd.concat([T, F.add_prefix("fail_")], axis=1); T["passed"] = DEC.passed(T); T["remote_evaluated"] = ev
    if BLIND:                                     # candidates' rows only in the sealed store (their positions and verdicts would unblind)
        T[T.population == "candidate"].to_parquet(f"{SEALED}/{pair}_candidates.parquet")
        T[T.population != "candidate"].to_parquet(ob + "_flags.parquet")
    else:
        T.to_parquet(ob + "_flags.parquet")
    flow = DEC.cutflow(T)
    n_dir = int(T[T.population == "chance"]["shift"].max() + 1) if (T.population == "chance").any() else 0
    G_half = SEL.efficiency_table(J, T[T.population == "injection"].set_index("inj").passed, ob) if "injection" in flow else np.nan
    r = dict(pair=pair, epoch=o["epoch"], candidates=flow.get("candidate", {}).get("rows", 0), final=flow.get("candidate", {}).get("passed", 0),
             chance_final=round(flow.get("chance", {}).get("passed", 0) / max(n_dir, 1), 3), chance_dirs=n_dir, injected=int(len(J)),
             chance_by_radius={int(A.background.chance_shift_arcsec): round(flow.get("chance", {}).get("passed", 0) / max(n_dir, 1), 3),
                               **{int(p_.split("_")[1]): round(flow[p_]["passed"] / A.background.shifts, 3) for p_ in flow if p_.startswith("chance_")}},
             injections_passed=flow.get("injection", {}).get("passed", 0), G_eff_half=G_half, cutflow=flow, source=SEL.CFG["source"],
             vetoes_complete=SEL.CFG["source"] == "primary", offline=SEL.CFG["offline"], rule=SEL.CFG["rule"],
             fidelity_source=fid_src,
             asteroid_ephemeris=SEL.asteroid_veto(C["pair"], np.zeros(0), np.zeros(0)) is not None,   # False: the asteroid cut is off for this pair
             scanner_marks=os.path.basename(SEL.CFG["scanner_marks"] or "none"),
             cache_used=sorted(os.path.basename(u) for u in SEL.USED))
    if BLIND:                                     # per-pair candidate counts are sealed; the open record keeps chance and injections
        json.dump(dict(pair=pair, candidates=r["candidates"], final=r["final"], candidate_cutflow=flow.get("candidate")), open(f"{SEALED}/{pair}_select.json", "w"), indent=1, default=float)
        r = DEC.open_record(r)
    json.dump(r, open(ob + "_select.json", "w"), indent=1, default=float); return r


if __name__ == "__main__":
    # --jobs pairs at once in spawned processes (the time is archive round trips, not computation; the answers do not
    # depend on the order or the concurrency)
    import multiprocessing as mp
    from concurrent.futures import ProcessPoolExecutor
    FILES = sorted(f for f in glob.glob(f"{PD}/*.json") if "_" not in os.path.basename(f) and (not ONLY or os.path.basename(f)[:-5] in ONLY))   # PAIR.json
    jobs = ARGS.jobs; rows = []
    from contextlib import nullcontext
    with ProcessPoolExecutor(jobs, mp_context=mp.get_context("spawn")) if jobs > 1 else nullcontext() as ex:
        for r in (ex.map(one, FILES) if ex else map(one, FILES)):
            rows.append(r)
            print(r["pair"], r.get("refused") or (("" if BLIND else f"candidates {r['candidates']} final {r['final']} ") + f"chance {r['chance_final']} injections {r['injections_passed']}/{r['injected']}"), flush=True)
    pd.DataFrame([{k: v for k, v in r.items() if k not in ("cutflow", "cache_used")} for r in rows]).to_csv(f"{OD}/select_summary.csv", index=False)
