"""The results report's numbers (notes/requirements.md U6), from the outputs after the unlock: the verdict, the two
tests, the fitted suppression and glint density, the upper limit as a share and as a rate, Bruehl et al.'s test as
published, the nulls, and the exploratory secondary from the seal (the anti-solar frame). Writes OUT.json with each number and the phrase the report states it in; the report's test
checks every phrase against the report's text, as P1 does for the pilot.
python -E -s -m heidelberg.producers.report_numbers VERDICT.json RATE.json OUT.json"""
import json
import os
import sys

V, R = (json.load(open(f)) for f in sys.argv[1:3]); OUT = sys.argv[3]
P = V["primary"]
limit = f"above {P['upper_limit_share']:.2f}" if P["upper_limit_above_grid"] else f"{P['upper_limit_share']:.2f}"
out = dict(verdict=V["verdict"], reasons=V["reasons"], likelihood_p=P["likelihood_p"], check_p=P["check_p"], kappa_hat=P["kappa_hat"],
           eta_hat=P["eta_hat"], share_95=P["upper_limit_share"], rate_95=R["rate_95"], bruehl_p=P["bruehl_p"], nulls=V["nulls"])
out["phrases"] = [
    f"The result is {V['verdict']}",
    f"the likelihood's p-value is {P['likelihood_p']:.3f} and the check's {P['check_p']:.3f}",
    f"the fitted suppression in shadow is kappa = {P['kappa_hat']:.2f}",
    f"the 95% upper limit on the glint share is {limit}",
    f"{R['rate_95']:.2g} glints per square degree per hour",
    f"Bruehl et al.'s test as published gives p = {P['bruehl_p']:.3f}",
] + [f"the {k.split('_')[0]} null test {'passes' if v['null_pass'] else 'fails'}" for k, v in V["nulls"].items()]
S = V.get("secondary_exploratory", {}); out["secondaries"] = S
out["phrases"] += [f"in the {k} (exploratory), the likelihood's p-value is {s['likelihood_p']:.3f} and the check's {s['check_p']:.3f}" for k, s in S.items()]
os.makedirs(os.path.dirname(os.path.abspath(OUT)), exist_ok=True)
json.dump(out, open(OUT, "w"), indent=1, default=float)
for p in out["phrases"]:
    print(p)
