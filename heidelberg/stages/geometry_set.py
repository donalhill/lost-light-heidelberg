"""The calibration's geometry: per pair of the record set (examined pairs left out, in the record set's order), the
shadow share on the 0.2 deg grid over the box from the record geometry (stages/record_geometry.py), and zero for the
sunlit pairs, which the shadow does not touch.
python -E -s -m heidelberg.stages.geometry_set RECORD_SET.csv GEOMETRY.csv GEOMETRY_s.npy OUT.npy"""
import sys

import numpy as np
import pandas as pd

RS, GC, GS, OUT = sys.argv[1:5]
R = pd.read_csv(RS); R = R[~R.examined]; G = pd.read_csv(GC); S = np.load(GS); idx = {a: i for i, a in enumerate(G.a)}
M = np.vstack([S[idx[a]] if role == "shadow" else np.zeros(S.shape[1]) for a, role in zip(R.a, R.role)])
assert (R.role == "shadow").sum() == int((M.max(1) > 0).sum()), "a shadow pair the geometry does not touch"
np.save(OUT, M); print(OUT, M.shape, f"{int((R.role == 'shadow').sum())} shadow, {int((R.role == 'sunlit').sum())} sunlit")
