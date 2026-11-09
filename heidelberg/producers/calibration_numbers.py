"""The calibration's numbers (requirement J3) as one JSON file for the draft's numbers (tools/numbers.py), from the text
calibration.sh prints: the sup critical value and its bootstrap standard error, the null's toys per share, the power
table (kappa 1 and 0.5; with refused pairs; with the times in error), the claim rule's and the within-pair test's false
rejection with the map correct, uniform, and at half the gradient, the likelihood's robustness cases, and the belt's
expected limit and coverage (also with the chance estimates scaled and with the times in error).
python -E -s -m heidelberg.producers.calibration_numbers CALIBRATION.txt OUT.json"""
import json
import re
import sys

import numpy as np

from heidelberg.config import A

TXT, OUT = sys.argv[1:3]
log = open(TXT).read()
blk = lambda name: log[log.index(name):]
unsup = lambda b, f: re.search(rf"glint share {f}, no shadow suppression: likelihood rejected [\d.]+, check rejected ([\d.]+), both ([\d.]+)", b).groups()
both = lambda b, f: float(unsup(b, f)[1])
size = lambda name: dict(no_glints=both(blk(name), "0.00"), glints_0p10=both(blk(name), "0.10"), glints_0p50=both(blk(name), "0.50"),
                         check_no_glints=float(unsup(blk(name), "0.00")[0]), check_glints_0p10=float(unsup(blk(name), "0.10")[0]),
                         check_glints_0p50=float(unsup(blk(name), "0.50")[0]))
table = lambda pre: {f"{float(f):.2f}": dict(likelihood=float(a), check=float(c), claim=float(b))
                     for f, a, c, b in re.findall(rf"^{pre}share ([\d.]+), likelihood ([\d.]+), check ([\d.]+), claim ([\d.]+)", log, re.M)}
limits = lambda pre: {f"{float(f):.2f}": dict(median=float(m), **({} if c is None or c == "" else dict(coverage=float(c))))
                      for f, m, c in re.findall(rf"^{pre}, share ([\d.]+): limit median ([\d.]+)(?:, covers in ([\d.]+))?", log, re.M)}
key = lambda x: f"{x:g}".replace(".", "p")
crit = float(re.search(r"sup critical value \(largest 95th percentile over glint shares\): ([\d.]+)", log).group(1))
q95 = {f: (float(q), float(se)) for f, q, se in re.findall(r"null_q_f([\d.]+)\.npy q95 ([\d.]+) bootstrap se ([\d.]+)", log)}
share_crit = max(q95, key=lambda f: q95[f][0])
power = {f"{float(f):.2f}_kappa{float(k):g}": dict(likelihood=float(a), check=float(c), claim=float(b))
         for f, k, a, c, b in re.findall(r"glint share ([\d.]+), kappa ([\d.]+): likelihood ([\d.]+), check ([\d.]+), both \(the claim\) ([\d.]+)", log)}
q50, q16, q84 = (float(x) for x in re.search(r"true share 0.00: belt upper limit median ([\d.]+) \(16-84%: ([\d.]+)-([\d.]+)\)", log).groups())
cover = {f"{float(f):.2f}": float(c) for f, c in re.findall(r"true share ([\d.]+): belt upper limit .*covers in ([\d.]+)", log)}
robust = {k.strip(): float(v) for k, v in re.findall(r"^(\S.*?)\s+rejected ([\d.]+) \(\d+ of \d+; nominal 0.05\)$", log, re.M)}
out = dict(tuning_background=float(re.search(r"tuning background per pair: ([\d.]+)", log).group(1)),
           critical_value=crit, critical_value_se=q95[share_crit][1], critical_share=float(share_crit),
           null_toys_per_share=int(re.search(r"glint share 0.1: .*\((\d+) toys\)", log).group(1)),
           null_shares=sorted(float(f) for f in q95), power=power,
           power_mc_se_max=max(float(np.sqrt(v * (1 - v) / A.calibration.power_toys)) for T in (power, table(r"pairs refused at [\d.]+: "), table(r"timing blur \(\d+ sets\): ")) for e in T.values() for v in e.values()),
           size_map_correct=size("claim rule, background map correct"), size_map_uniform=size("claim rule, background uniform"),
           size_map_half=size("claim rule, background map half"),
           refused_share=float(re.search(r"^pairs refused at ([\d.]+):", log, re.M).group(1)), power_refused=table(r"pairs refused at [\d.]+: "),
           blur_sets=int(re.search(r"^timing blur \((\d+) sets\):", log, re.M).group(1)), power_blur=table(r"timing blur \(\d+ sets\): "),
           limit_chance_scale={key(s): limits(f"chance scale {s:.2f}") for s in A.calibration.cal_chance_scales}, limit_blur=limits("timing blur"),
           robustness=robust, expected_limit=dict(median=q50, lo16=q16, hi84=q84), coverage=cover,
           belt_toys=int(re.search(r"true share 0.00: belt upper limit .*\((\d+) toys\)", log).group(1)))
json.dump(out, open(OUT, "w"), indent=1); print(f"{OUT}: critical value {crit} (share {share_crit}), expected limit {q50}")
