"""The sample's numbers in the pre-registration (notes/traceability.md G2), from GAVO's metadata alone: per record pair
(examined pairs left out) the plates' start times, centres, and exposures; the shadow's coverage of the test area
(record geometry); the pairs removed as not eligible (requirement S7: the touched pairs outside the record set, the
tuning and examined pairs aside, and the sunlit pairs short of the 113 matched) and the plates whose times moved from
GAVO's dateObs (S6); the medians and standardised differences of the covariates by role (record set); the test area on
the tuning pairs (the dry run's coverage map); the plates on which GAVO's dateObs is 12 h early (all with exposures
running past 0 h UT), and the logbook plates (producers/logbook.py). LSW.csv is the retimed plate table
(stages/plate_times.py), LSW_RAW.csv GAVO's. Writes OUT.json with each number and the phrase the draft states it in;
requirement S5 checks every phrase against the draft.
python -E -s -m heidelberg.producers.sample_numbers LSW.csv RECORD_SET.csv GEOMETRY.csv TUNING_COVERAGE.csv ALL_PAIRS.csv TUNING_PAIRS.csv LSW_RAW.csv OUT.json"""
import json
import os
import sys

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord

from heidelberg.stages.role_report import std_diff
from heidelberg.producers.logbook import LOGBOOK

from heidelberg.config import A

LSW, RS, GC, COV, ALL, TUN, RAW, OUT = sys.argv[1:9]; S_ = A.sample; TOL_MIN = A.times.dateobs_tol_min; PER = S_.sunlit_per_shadow; Y49 = S_.poss_first_year
H12 = json.load(open(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "heidelberg", "facts.json")))["gavo_offset_h"]["value"]
L = pd.read_csv(LSW).drop_duplicates("plateID").set_index("plateID"); R0 = pd.read_csv(RS); R = R0[~R0.examined].reset_index(drop=True)
G = pd.read_csv(GC).set_index("a"); C = pd.read_csv(COV)
Pa, Pb = L.loc[R.a], L.loc[R.b]
dt = np.abs((pd.to_datetime(Pa.startTime.values) - pd.to_datetime(Pb.startTime.values)).total_seconds())
sep = SkyCoord(Pa.centerAlpha.values, Pa.centerDelta.values, unit="deg").separation(SkyCoord(Pb.centerAlpha.values, Pb.centerDelta.values, unit="deg")).arcmin
eq = float(np.mean(Pa.exposure.values == Pb.exposure.values))
tun = set(pd.read_csv(TUN).a); touched = [a for a in G.index if a not in tun and a not in set(R0.a[R0.examined])]
Lc = L[L.time_source == "catalogue"]; Raw = pd.read_csv(RAW).drop_duplicates("plateID").set_index("plateID").loc[Lc.index]
twelve = ((Raw.dateObs - Lc.dateObs) * 1440 + H12 * 60).abs() <= TOL_MIN   # GAVO's dateObs 12 h before the catalogue's
t0 = pd.to_datetime(Lc.startTime); past0 = (t0 + pd.to_timedelta(Lc.span_s, unit="s")).dt.date > t0.dt.date
assert past0[twelve].all(), "a plate 12 h early in GAVO whose exposure does not run past 0 h UT"
years = sorted({d[:4] for _, d, *_ in LOGBOOK})
n12_early = int((t0[twelve].dt.year < Y49).sum())
removed_shadow = len(touched) - int((R.role == "shadow").sum())
drawn_sunlit = round(PER * sum(a not in tun for a in G.index))   # as stages/record_set.py draws: the examined pair included
removed_sunlit = drawn_sunlit - int((R.role == "sunlit").sum())
sh = R.role == "shadow"; Gs = G.loc[R.a[sh]]
# the central share's percentiles are rounded: 50 x (1 - 0.9) is 4.999999999999999 in floating point
out = dict(pairs=len(R), shadow=int(sh.sum()), sunlit=int((~sh).sum()), dt_median=float(np.median(dt)), dt_95=float(np.percentile(dt, 100 * S_.dt_share)),
           sep_median=float(np.median(sep)), sep_90=[float(np.percentile(sep, round(50 * (1 - S_.sep_share), 10))), float(np.percentile(sep, round(50 * (1 + S_.sep_share), 10)))], exposure_equal=eq,
           edge_crossing=int((Gs.partial > S_.edge_area).sum()), mostly=int((Gs.full >= S_.mostly_share).sum()), wholly=int((Gs.full >= 1.0 - 1e-12).sum()),
           neither=int(len(Gs) - (Gs.partial > S_.edge_area).sum() - (Gs.full >= S_.mostly_share).sum()),   # the two classes cannot overlap
           removed_shadow=removed_shadow, removed_sunlit=removed_sunlit, sunlit_drawn=drawn_sunlit, plates_eligible=int(L.eligible.sum()), plates=len(L),
           gavo_12h=int(twelve.sum()), gavo_12h_before_1949=n12_early, logbook_plates=len(LOGBOOK), logbook_years=years,
           antisun_median=[float(R.antisun_deg[sh].median()), float(R.antisun_deg[~sh].median())],
           ecl_lat_median=[float(R.ecl_lat[sh].median()), float(R.ecl_lat[~sh].median())],
           ecl_lat_std_diff=std_diff(R.ecl_lat[sh], R.ecl_lat[~sh]), area=[float(C.area_deg2.mean()), float(C.area_deg2.std())])
AP = pd.read_csv(ALL); late = AP[AP.year >= Y49]; touched_late = int((G.year >= Y49).sum())
words = ["no", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine"]
out.update(pairs_all=len(AP), pairs_1949=len(late), touched_1949=touched_late, window_first_year=int(str(S_.first_date)[:4]),
           year_first=int(R.year.min()), year_last=int(R.year.max()),
           shadow_year_first=int(R.year[sh].min()), shadow_year_last=int(R.year[sh].max()),
           shadow_year_last_before_poss=int(R.year[sh & (R.year < Y49)].max()))
last = pd.Timestamp("1858-11-17") + pd.to_timedelta(float(R.mjd[sh].max()), unit="D")   # MJD to UT
SPLIT = 1935                                                      # where the archive's twin exposures fall away (pairs per year)
out.update(archive_split_year=SPLIT, pairs_before_split=int((AP.year < SPLIT).sum()), pairs_from_split=int((AP.year >= SPLIT).sum()))
out.update(shadow_last_date=f"{last.day} {last:%B %Y}", shadow_last_days_before_end=int((pd.Timestamp(str(S_.end_before)) - last.normalize()).days))
sign = lambda x, d: ("−" if x < 0 else "") + f"{abs(x):.{d}f}"
out["phrases"] = [
    f"Of the {len(late)} Bruce twin pairs exposed between {Y49} and October 1957, the years of POSS-I, Earth's shadow touches the central field of only {words[touched_late]}",
    f"There are {len(AP):,} such pairs",
    f"start a median of {out['dt_median']:.0f} seconds apart (95% within {out['dt_95']:.0f} seconds)",
    f"their centres lie a median of {out['sep_median']:.1f} arcminutes apart",
    f"{out['sep_90'][0]:.1f} to {out['sep_90'][1]:.1f} arcminutes for 90% of pairs",
    f"their exposures are equal on {100 * eq:.0f}% of pairs",
    f"The record set is the {out['shadow']} pairs whose test area Earth's shadow touches, together with {out['sunlit']} sunlit pairs",
    f"Of the {out['shadow']} shadow pairs, {out['edge_crossing']} have the shadow's edge (where a point is in shadow for between 5% and 95% of the exposure) over more than a tenth of the test area, {out['mostly']} have at least 95% of it in full shadow, and {out['wholly']} lie wholly in shadow",
    f"The test areas are nearly equal, {out['area'][0]:.1f} ± {out['area'][1]:.1f} square degrees on the tuning pairs",
    f"On {({12: 'twelve'} | dict(enumerate(words))).get(len(LOGBOOK), len(LOGBOOK))} plates from {', '.join(years[:-1])}, and {years[-1]} the catalogue's start and end agree with the scanned observation logbooks",
    f"GAVO's published mid-exposure is 12 hours early on {twelve.sum():,} plates whose exposure ran past 0 h UT ({words[n12_early] if n12_early < 10 else n12_early} of them before {Y49})",
    f"removes {words[removed_shadow] if removed_shadow < 10 else removed_shadow} shadow pairs and {words[removed_sunlit] if removed_sunlit < 10 else removed_sunlit} sunlit pairs",
    f"distance from the anti-solar point (the latter's medians {out['antisun_median'][0]:.1f} and {out['antisun_median'][1]:.1f} degrees for shadow and sunlit pairs)",
    f"not in ecliptic latitude (standardised difference of the signed latitude {sign(out['ecl_lat_std_diff'], 2)})",
]
json.dump(out, open(OUT, "w"), indent=1, default=float)
print(json.dumps({k: v for k, v in out.items() if k != "phrases"}, default=float))
for p in out["phrases"]:
    print(" ", p)
