"""Data quality by role before unblinding (notes/requirements.md C2, C3). For each pair of the record set, quantities that
carry no candidate result: whether both solutions were accepted, whether the pair passed the selection (and why not),
depth (G50, the mean over the plates of the G at which half of the Gaia stars are detected), aperture noise, the
efficiency E_i (lib/efficiency.py), the G at half the peak efficiency, the chance estimate, the recurrence coverage of
its null locations, and from GAVO's metadata whether both plates' times are self-consistent (dateObs within 10 minutes of
start + exposure / 2) and the twins' common exposure. Shadow and sunlit pairs are compared by the standardised
difference (difference of the means over the root of the mean of the variances). The role rule: a standardised
difference beyond 0.5 in acceptance, selection, depth, aperture noise, or E_i means no claim. Also printed: by role, the
medians of ecliptic latitude and of the distance from the anti-solar point, from the record set.
python -E -s -m heidelberg.stages.role_report RECORD_SET.csv PRODUCTS_DIR SELECT_DIR RECURRENCE_DIR|- LSW.csv|- OUT.csv"""
import json
import os
import sys

import numpy as np
import pandas as pd

from heidelberg import config
from heidelberg.lib.efficiency import pair_efficiency

RULE = ("accepted", "selected", "g50", "apnoise", "eff")
TCON = config.A.roles.time_consistent_min            # min: GAVO's mid-exposure against start + exposure / 2
REPORTED = ("accepted", "selected", "g50", "apnoise", "eff", "G_eff_half", "detections", "chance", "null_informative", "time_consistent", "overlap")


def std_diff(a, b):
    a, b = pd.to_numeric(a, errors="coerce").dropna(), pd.to_numeric(b, errors="coerce").dropna()
    sd = np.sqrt(0.5 * (a.var() + b.var())) if len(a) > 1 and len(b) > 1 else 0.0
    return float((a.mean() - b.mean()) / sd) if sd > 0 else 0.0


def rule(D):
    """The quantities of RULE whose standardised difference between roles exceeds 0.5; a quantity missing from the
    report fails (a broken report cannot pass the rule)."""
    S, U = D[D.role == "shadow"], D[D.role == "sunlit"]
    return [v if v in D else f"{v} (missing)" for v in RULE if v not in D or abs(std_diff(S[v], U[v])) > config.A.roles.role_rule]


def table(R, PD, SD, RD, LSW):
    TM = None
    if LSW:
        from astropy.time import Time
        Lw = pd.read_csv(LSW).dropna(subset=["startTime"])
        Lw["t_off_min"] = (Lw.dateObs - (Time(list(Lw.startTime), format="isot", scale="utc").mjd + Lw.exposure / 172800)) * 1440
        TM = Lw.set_index("plateID")
    REC = pd.read_csv(f"{RD}/recurrence_locations.csv") if RD and os.path.exists(f"{RD}/recurrence_locations.csv") else None
    rows = []
    for _, p in R.iterrows():
        r = dict(pair=p.a, role=p.role, s_mean=p.s_mean, measured=os.path.exists(f"{PD}/{p.a}.json"))
        if TM is not None and p.a in TM.index and p.b in TM.index:
            A, B = TM.loc[p.a], TM.loc[p.b]
            lo = max(A.dateObs - A.exposure / 172800, B.dateObs - B.exposure / 172800); hi = min(A.dateObs + A.exposure / 172800, B.dateObs + B.exposure / 172800)
            r.update(time_consistent=bool(abs(A.t_off_min) <= TCON and abs(B.t_off_min) <= TCON), overlap=float((hi - lo) * 86400 / min(A.exposure, B.exposure)))
        if r["measured"]:
            o = json.load(open(f"{PD}/{p.a}.json")); sv = o["solve"]
            g50 = [x.get("g50_astrometric", np.nan) for x in sv]
            r.update(accepted=all(x.get("accepted") for x in sv), g50=float(np.nanmean([g if np.isscalar(g) else np.nan for g in g50])),
                     apnoise=float(np.mean(o["aperture_noise"])), detections=float(np.mean(o["detections"])))
        if os.path.exists(f"{SD}/{p.a}_select.json"):
            s = json.load(open(f"{SD}/{p.a}_select.json"))
            r.update(refused=str(s.get("refused", ""))[:40], selected="refused" not in s, chance=s.get("chance_final", np.nan), G_eff_half=s.get("G_eff_half", np.nan))
            fl = f"{SD}/{p.a}_flags.parquet"
            if "refused" not in s and os.path.exists(fl):
                r["eff"] = pair_efficiency(pd.read_parquet(f"{PD}/{p.a}_inj.parquet"), pd.read_parquet(fl, columns=["population", "passed", "inj"]))
        if REC is not None:
            d = REC[(REC.pair == p.a) & (REC.kind == "blank")]
            r.update(null_locations=len(d), null_informative=float(d.n_inf.mean()) if len(d) else np.nan)
        rows.append(r)
    return pd.DataFrame(rows)


if __name__ == "__main__":
    RS, PD, SD, RD, LSW, OUT = sys.argv[1:7]; RD = None if RD == "-" else RD; LSW = None if LSW == "-" else LSW
    R = pd.read_csv(RS); R = R[~R.examined] if "examined" in R else R
    D = table(R, PD, SD, RD, LSW); D.to_csv(OUT, index=False)
    S, U = D[D.role == "shadow"], D[D.role == "sunlit"]
    print(f"{len(D)} pairs ({len(S)} shadow, {len(U)} sunlit); measured {int(D.measured.sum())}")
    print(f"{'quantity':18s} {'shadow':>10s} {'sunlit':>10s} {'std diff':>9s}")
    for v in REPORTED:
        if v in D and D[v].notna().any():
            print(f"{v:18s} {pd.to_numeric(S[v], errors='coerce').mean():10.3f} {pd.to_numeric(U[v], errors='coerce').mean():10.3f} {std_diff(S[v], U[v]):9.2f}")
    for v in ("ecl_lat", "antisun_deg"):
        if v in R:
            print(f"{v + ' (median)':18s} {R[R.role == 'shadow'][v].median():10.3f} {R[R.role == 'sunlit'][v].median():10.3f}")
    bad = rule(D)
    print(f"role rule (|std diff| <= {config.A.roles.role_rule} on " + ", ".join(RULE) + "): " + ("PASS" if not bad else "FAIL on " + ", ".join(bad)))
    if "refused" in D:
        print("refusals by role:", D[D.refused.astype(str) != ""].groupby("role").refused.value_counts().to_dict())
