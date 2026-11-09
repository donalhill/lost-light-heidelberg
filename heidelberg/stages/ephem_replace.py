"""Replace some plates' asteroid ephemerides in an ephemeris folder with ones recomputed at new times (plates whose
times were corrected after their ephemerides were computed), or add plates new to the folder. Every PREFIX_part*.parquet
in EPH_DIR loses those plates' rows; NEW_DIR's parts are added as EPH_DIR/PREFIX_retimed_part*.parquet; plates.csv takes
the new times and exposures, a plate new to it appended. The selection reads every parquet in the folder, so no stale
row may remain.
python -E -s -m heidelberg.stages.ephem_replace EPH_DIR NEW_DIR TIMES.csv PREFIX"""
import glob
import os
import sys

import pandas as pd

EPH, NEW, TIMES, PREFIX = sys.argv[1:5]
T = pd.read_csv(TIMES); plates = set(T.plate); dropped = 0
for f in sorted(glob.glob(f"{EPH}/*_part*.parquet")):
    A = pd.read_parquet(f)
    if "plate" in A and A.plate.isin(plates).any():
        dropped += int(A.plate.isin(plates).sum()); A[~A.plate.isin(plates)].to_parquet(f)
added = 0
for f in sorted(glob.glob(f"{NEW}/asteroids_part*.parquet")):
    A = pd.read_parquet(f)
    if "plate" not in A:                                          # a batch with no asteroid near these fields
        continue
    assert set(A.plate) <= plates, f"{f} holds other plates"
    added += len(A); A.to_parquet(f"{EPH}/{PREFIX}_retimed_{os.path.basename(f).split('_', 1)[1]}")
P = pd.read_csv(f"{EPH}/plates.csv").set_index("plate"); t = T.set_index("plate"); old = t.index.intersection(P.index, sort=False)
P.loc[old, "ut"] = t.ut[old]; P.loc[old, "exp"] = t.exp[old]
new = t.loc[t.index.difference(P.index, sort=False), ["ut", "exp"]]; P = pd.concat([P, new])
P.rename_axis("plate").reset_index().to_csv(f"{EPH}/plates.csv", index=False)
print(f"{len(plates)} plates ({len(new)} new to the folder): {dropped} stale rows removed, {added} rows added; plates.csv updated")
