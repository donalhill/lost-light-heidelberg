"""Whether a rebuilt product equals the one in use (requirement J5): CSV tables by value (the same columns and rows;
every non-float entry equal; floats equal within 1e-9, absolute or relative, since another platform's maths library
may differ in the last digits), .npy arrays likewise, other files byte for byte. The largest float difference is
printed. With --drop-examined the second CSV is the first without its examined rows and column (the record set against
heidelberg/record_pairs.csv).
python -E -s tools/same.py NEW REFERENCE [--drop-examined]"""
import filecmp
import sys

import numpy as np
import pandas as pd

RTOL, ATOL = 1e-9, 1e-9              # |a - b| <= ATOL + RTOL |b|: shares near 0 differ by ~1e-10 across platforms


def close(a, b):
    """Whether two arrays are equal (floats within ATOL + RTOL |b|, NaN equal to NaN), and their largest float difference."""
    a, b = np.asarray(a), np.asarray(b)
    if a.shape != b.shape:
        return False, np.inf
    if a.dtype.kind != "f" and b.dtype.kind != "f":
        return bool(pd.Series(a.ravel()).equals(pd.Series(b.ravel()))), 0.0      # missing entries equal
    if a.dtype.kind not in "fiu" or b.dtype.kind not in "fiu":
        return False, np.inf
    a, b = a.astype(float), b.astype(float)
    d = float(np.nanmax(np.abs(a - b), initial=0.0)) if a.size else 0.0
    return bool(np.allclose(a, b, rtol=RTOL, atol=ATOL, equal_nan=True)), d


def same(new, ref, drop_examined=False):
    if new.endswith(".csv"):
        A, B = pd.read_csv(new), pd.read_csv(ref)
        if drop_examined:
            A = A[~A.examined].drop(columns="examined").reset_index(drop=True)
        if list(A.columns) != list(B.columns) or len(A) != len(B):
            return False, np.inf
        res = [close(A[c].to_numpy(), B[c].to_numpy()) for c in A.columns]
        return all(r[0] for r in res), max((r[1] for r in res), default=0.0)
    if new.endswith(".npy"):
        return close(np.load(new), np.load(ref))
    return filecmp.cmp(new, ref, shallow=False), 0.0


if __name__ == "__main__":
    new, ref = sys.argv[1:3]
    ok, d = same(new, ref, "--drop-examined" in sys.argv)
    print(("same: " if ok else "DIFFERENT: ") + f"{new} and {ref}" + (f" (largest float difference {d:.3g})" if d else ""))
    sys.exit(0 if ok else 1)
