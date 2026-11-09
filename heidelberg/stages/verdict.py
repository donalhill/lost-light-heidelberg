"""The verdict after the unlock (notes/requirements.md U4), the decision rule of the pre-registration applied by code
frozen in the registration, one of three outcomes: a claim of a deficit in Earth's shadow if the likelihood and the check
both reject H0 at the 5% level (the primary selection, topocentric frame, area controls), both null tests pass, and the
comparison of data quality between roles passes; null if neither test rejects and every check passes (the upper limit
is then the result); inconclusive in every other case (one test rejecting and not the other, a null test failing, or
the comparison of data quality failing), the upper limit reported beside the reason it is in doubt. Every reason is
listed. The secondary analysis (the anti-solar frame) is reported beside the result, labelled exploratory, and
decides nothing.
python -E -s -m heidelberg.stages.verdict SHADOW_DIR ROLE_REPORT.csv OUT.json
  SHADOW_DIR: the unsealed shadow outputs (sealed_candidate_*.json) with the open nulls (null_*.json)."""
import json
import os
import sys

import pandas as pd

from heidelberg.stages.role_report import rule as role_rule

from heidelberg.config import A

ALPHA = A.statistics.alpha
NULLS = ("null_chance_area_topo.json", "null_injection_area_topo.json")
PRIMARY = "sealed_candidate_area_topo.json"
SECONDARY = {"anti-solar frame": "sealed_candidate_area_geo.json"}


OUTCOMES = {"claim": "a deficit in Earth's shadow is claimed", "null": "null: no deficit in Earth's shadow; the upper limit is the result",
            "inconclusive": "inconclusive: no claim; the upper limit rests on a background model that the reasons call into question"}


def decide(primary, nulls, role_failed):
    """primary: the primary shadow output; nulls: {name: null output}; role_failed: the quantities failing the role
    rule. Returns (outcome, reasons): claim, null, or inconclusive, and every condition that fails."""
    lik, chk = primary["likelihood"]["p"] < ALPHA, primary["p"] < ALPHA; reasons = []
    if not lik:
        reasons.append(f"the likelihood does not reject (p = {primary['likelihood']['p']:.4f})")
    if not chk:
        reasons.append(f"the check does not reject (p = {primary['p']:.4f})")
    failed = [f"the null test {name} fails (likelihood p = {n['likelihood']['p']:.4f}, check p = {n['p']:.4f} and {n['p_excess']:.4f})"
              for name, n in nulls.items() if not n["null_pass"]]
    if role_failed:
        failed.append("the comparison of data quality between roles fails on " + ", ".join(role_failed))
    if lik and chk and not failed:
        return "claim", []
    if not lik and not chk and not failed:
        return "null", reasons
    return "inconclusive", reasons + failed


def summary(o):
    L = o["likelihood"]; U = L.get("upper_limit", {})
    return dict(likelihood_p=L["p"], check_p=o["p"], kappa_hat=L["kappa_hat"], eta_hat=L.get("eta_hat"), scale_hat=L.get("scale_hat"),
                upper_limit_share=U.get("share_95"), upper_limit_above_grid=U.get("above_grid"), upper_limit_below_grid=U.get("below_grid"), rate_per_pair=U.get("rate_per_pair"),
                finals=o["N"], bruehl_p=o["bruehl_replication"]["p"])


if __name__ == "__main__":
    SD, RR, OUT = sys.argv[1:4]
    primary = json.load(open(f"{SD}/{PRIMARY}")); nulls = {n[5:-5]: json.load(open(f"{SD}/{n}")) for n in NULLS}
    failed = role_rule(pd.read_csv(RR))
    outcome, reasons = decide(primary, nulls, failed)
    out = dict(outcome=outcome, claimed=outcome == "claim", verdict=OUTCOMES[outcome], reasons=reasons,
               primary=summary(primary), nulls={k: dict(null_pass=v["null_pass"]) for k, v in nulls.items()}, role_rule_failed=failed,
               secondary_exploratory={k: summary(json.load(open(f"{SD}/{f}"))) for k, f in SECONDARY.items() if os.path.exists(f"{SD}/{f}")})
    json.dump(out, open(OUT, "w"), indent=1, default=float)
    print(out["verdict"] + ("" if outcome == "claim" else ": " + "; ".join(reasons)))
