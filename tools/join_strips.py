"""The VSX strips fetched from VizieR joined into one table in the layout of the other VizieR snapshots (a header line,
the units and dashes rows, then tab-separated rows with _RAJ2000 and _DEJ2000), each object once (a strip's edge
belongs to both of its strips), and its row count printed.
python tools/join_strips.py STRIP_DIR OUT.tsv"""
import glob
import sys

import pandas as pd

parts = []
for f in sorted(glob.glob(f"{sys.argv[1]}/*.tsv")):
    T = pd.read_csv(f, sep="\t", comment="#", skip_blank_lines=True, dtype=str)
    head, T = T.iloc[:2], T.iloc[2:]
    if len(T) >= 1_000_000 - 1:
        sys.exit(f"{f}: {len(T)} rows, at VizieR's row cap: split the strip")
    parts.append(T)
V = pd.concat(parts).drop_duplicates("OID")
pd.concat([head, V]).to_csv(sys.argv[2], sep="\t", index=False)
print(f"{len(V)} VSX objects from {len(parts)} strips")
