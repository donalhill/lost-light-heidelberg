"""The transient sample's checksum (notes/requirements.md R1, R2): the SHA-256 of a canonical table of the passing
candidates of a blind selection (SELECT_DIR/sealed/*_candidates.parquet): per transient its pair and the two
detections' positions to 1e-7 deg, sorted, as fixed-format text. Nothing about the sample is printed but the hash (its
size would tell the number of candidates).
python -E -s -m heidelberg.stages.sample_checksum SELECT_DIR [REGISTERED_SHA256]
  with REGISTERED_SHA256 the run exits 1 unless the sample's checksum equals it (part two's check)."""
import glob
import hashlib
import os
import sys

import pandas as pd

COLS = ("ra", "dec", "ra_bo", "dec_bo")


def canonical(sd):
    rows = []
    for f in sorted(glob.glob(f"{sd}/sealed/*_candidates.parquet")):
        T = pd.read_parquet(f); T = T[T.passed]
        pair = os.path.basename(f)[:-len("_candidates.parquet")]
        rows += [(pair,) + tuple(round(float(v), 7) for v in r) for r in T[list(COLS)].itertuples(index=False)]
    return "pair," + ",".join(COLS) + "\n" + "".join(f"{p}," + ",".join(f"{v:.7f}" for v in r) + "\n" for p, *r in sorted(rows))


def checksum(sd):
    return hashlib.sha256(canonical(sd).encode()).hexdigest()


if __name__ == "__main__":
    h = checksum(sys.argv[1]); print(h)
    if len(sys.argv) > 2:
        ok = h == sys.argv[2].strip(); print("the transient sample is as registered" if ok else "the transient sample DIFFERS from the registration"); sys.exit(0 if ok else 1)
