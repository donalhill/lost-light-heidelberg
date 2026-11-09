"""Known asteroids on each plate, computed ourselves (SkyBoT has been unreachable): MPCORB osculating elements
integrated back to the plate epochs with ASSIST (DE440 planets + the 16 massive asteroids, GR), seen from Heidelberg-Koenigstuhl (024) with
light time, at the start, middle, and end of each exposure. One batch of orbits per process; every batch steps through all
plate times in order, so each orbit is integrated once.
ONLY_DESIG (environment): a file of designations; only those orbits are integrated (for a plate whose time moved by
less than the margin kept around each field allows, the asteroids already found near it). The settings are in
config/analysis.json (section ephemerides) and config/compute.json, the physical constants in heidelberg/facts.json.
python -E -s -m heidelberg.stages.asteroid_ephem_024 PLATE_TIMES.csv DATA_DIR OUTDIR
-> OUTDIR/asteroids_part<k>.parquet: plate, number/designation, H, U, ra/dec at start/mid/end, V, r, delta"""
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from heidelberg import config

F, E_, CP = config.facts(), config.A.ephemerides, config.C.processes
K = F["gauss_k"]                        # Gauss, AU^1.5 / day
EPS = np.radians(F["obliquity_j2000_arcsec"] / 3600)   # J2000 obliquity, as MPC
C_AU_D = F["c_au_day"]                  # speed of light, AU / day
PH = config.A.physics                   # the site, as the shadow model's (Koenigstuhl; MPC 024)
LON, LAT, HGT = np.radians(PH.site_lon_deg), np.radians(PH.site_lat_deg), PH.site_height_m / 1000   # height in km
RE_AU = F["wgs84_a_km"] / F["au_km"]
FIELD = E_.ephem_field_deg              # deg from the plate centre, as SkyBoT cone


def unpack_epoch(s):
    """MPC packed date (e.g. K25BL) -> JD (TT)."""
    yr = {"I": 1800, "J": 1900, "K": 2000}[s[0]] + int(s[1:3])
    dig = lambda c: int(c) if c.isdigit() else ord(c) - 55
    return pd.Timestamp(year=yr, month=dig(s[3]), day=dig(s[4])).to_julian_date()


def read_mpcorb(path, hmax):
    cols = [(0, 7), (8, 13), (14, 19), (20, 25), (26, 35), (37, 46), (48, 57), (59, 68), (70, 79), (92, 103), (105, 106), (166, 194)]
    names = ["desig", "H", "G", "epoch", "M", "peri", "node", "incl", "e", "a", "U", "name"]
    import gzip
    with gzip.open(path, "rt", errors="replace") as f:
        skip = next(i for i, line in enumerate(f) if line.startswith("-----")) + 1
    D = pd.read_fwf(path, colspecs=cols, names=names, skiprows=skip, compression="infer", dtype=str, encoding_errors="replace")
    D = D[D.a.notna() & D.H.notna()]
    for c in ["H", "G", "M", "peri", "node", "incl", "e", "a"]:
        D[c] = pd.to_numeric(D[c], errors="coerce")
    D = D[(D.H <= hmax) & (D.e < 1) & D.a.notna()].copy()
    D["G"] = D.G.fillna(F["hg_default_g"])
    D["jd"] = D.epoch.map({e: unpack_epoch(e) for e in D.epoch.unique()})
    return D.reset_index(drop=True)


def elements_to_state(D):
    """Heliocentric ecliptic osculating elements -> heliocentric equatorial (ICRF) position and velocity, AU and AU/day."""
    a, e = D.a.values, D.e.values
    i, O, w, M = (np.radians(D[c].values) for c in ["incl", "node", "peri", "M"])
    E = M.copy()
    for _ in range(E_.ephem_kepler_iter):
        E -= (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
    n = K / a ** 1.5
    xp, yp = a * (np.cos(E) - e), a * np.sqrt(1 - e ** 2) * np.sin(E)
    Ed = n / (1 - e * np.cos(E))
    vxp, vyp = -a * np.sin(E) * Ed, a * np.sqrt(1 - e ** 2) * np.cos(E) * Ed
    cO, sO, cw, sw, ci, si = np.cos(O), np.sin(O), np.cos(w), np.sin(w), np.cos(i), np.sin(i)
    P = np.array([cO * cw - sO * sw * ci, sO * cw + cO * sw * ci, sw * si])
    Q = np.array([-cO * sw - sO * cw * ci, -sO * sw + cO * cw * ci, cw * si])
    r = P * xp + Q * yp; v = P * vxp + Q * vyp
    rot = np.array([[1, 0, 0], [0, np.cos(EPS), -np.sin(EPS)], [0, np.sin(EPS), np.cos(EPS)]])
    return (rot @ r).T, (rot @ v).T


def observer(jd_utc, earth):
    """Barycentric ICRF position of the site: the Earth's centre plus the site vector (equator of date; the 1950s precession
    shifts that 6,400 km vector by under 1%, a sub-arcsecond effect even for near-Earth objects at 0.05 AU)."""
    T = (jd_utc - 2451545.0) / 36525
    g0, g1, g2 = F["gmst_1982"]
    gmst = np.radians((g0 + g1 * (jd_utc - 2451545.0) + g2 * T ** 2) % 360)
    lst = gmst + LON
    f = 1 / F["wgs84_inv_flattening"]; C = 1 / np.sqrt(np.cos(LAT) ** 2 + (1 - f) ** 2 * np.sin(LAT) ** 2); Sf = (1 - f) ** 2 * C
    rc = (C + HGT / F["wgs84_a_km"]) * np.cos(LAT) * RE_AU; rs = (Sf + HGT / F["wgs84_a_km"]) * np.sin(LAT) * RE_AU
    return earth + np.array([rc * np.cos(lst), rc * np.sin(lst), rs])


def run_batch(args):
    k, D, times, data, out = args
    import assist
    import rebound
    for f in ("de440.bsp", "sb441-n16.bsp"):                 # ASSIST silently drops the asteroid forces without the second
        if not (Path(data) / f).exists():
            sys.exit(f"{data}/{f} missing")
    ephem = assist.Ephem(str(Path(data) / "de440.bsp"), str(Path(data) / "sb441-n16.bsp"))
    t0 = time.time()
    r, v = elements_to_state(D)
    rows = []
    # orbits sharing an epoch start together (MPCORB has a handful of epochs); integrate each group through every time
    states = np.zeros((len(D), 6, len(times)))
    for jd, idx in D.groupby("jd").indices.items():
        sim = rebound.Simulation()
        sun = ephem.get_particle("sun", jd - ephem.jd_ref)
        for j in idx:
            sim.add(x=r[j, 0] + sun.x, y=r[j, 1] + sun.y, z=r[j, 2] + sun.z, vx=v[j, 0] + sun.vx, vy=v[j, 1] + sun.vy, vz=v[j, 2] + sun.vz)
        sim.t = jd - ephem.jd_ref
        sim.ri_ias15.min_dt = E_.ephem_min_dt_d
        ex = assist.Extras(sim, ephem)                    # held: the forces stay attached to sim while it exists
        buf = np.zeros((len(idx), 3)); bv = np.zeros((len(idx), 3))
        for q, t in sorted(enumerate(times.tdb.values), key=lambda z: -z[1]):   # backwards in time from the epoch
            sim.integrate(t - ephem.jd_ref)
            sim.serialize_particle_data(xyz=buf, vxvyvz=bv)
            states[idx, :3, q] = buf; states[idx, 3:, q] = bv
    for q, tm in enumerate(times.itertuples()):
        tt = tm.tdb - ephem.jd_ref
        e = ephem.get_particle("earth", tt); s = ephem.get_particle("sun", tt)
        o = observer(tm.jd_utc, np.array([e.x, e.y, e.z]))
        p, pv = states[:, :3, q], states[:, 3:, q]
        d = p - o; tau = np.linalg.norm(d, axis=1) / C_AU_D
        d = d - pv * tau[:, None]                                        # light time (one iteration: tau < 0.1 d)
        dist = np.linalg.norm(d, axis=1)
        ra = np.degrees(np.arctan2(d[:, 1], d[:, 0])) % 360; dec = np.degrees(np.arcsin(d[:, 2] / dist))
        cd = np.sin(np.radians(dec)) * np.sin(np.radians(tm.dec_c)) + np.cos(np.radians(dec)) * np.cos(np.radians(tm.dec_c)) * np.cos(np.radians(ra - tm.ra_c))
        m = cd > np.cos(np.radians(FIELD + E_.ephem_margin_deg))
        if not m.any():
            continue
        hs = p[m] - np.array([s.x, s.y, s.z]); rh = np.linalg.norm(hs, axis=1)
        ph = np.arccos(np.clip(np.sum(hs * d[m], 1) / (rh * dist[m]), -1, 1))
        t2 = np.tan(ph / 2); G = D.G.values[m]
        A1, B1, A2, B2 = F["hg_phase"]
        phi = (1 - G) * np.exp(-A1 * t2 ** B1) + G * np.exp(-A2 * t2 ** B2)
        V = D.H.values[m] + 5 * np.log10(rh * dist[m]) - 2.5 * np.log10(phi)
        rows.append(pd.DataFrame({"plate": tm.plate, "when": tm.when, "desig": D.desig.values[m], "name": D.name.values[m],
                                  "H": D.H.values[m], "U": D.U.values[m], "ra": ra[m], "dec": dec[m], "V": V, "r": rh, "delta": dist[m]}))
    R = pd.concat(rows) if rows else pd.DataFrame()
    R.to_parquet(Path(out) / f"asteroids_part{k:04d}.parquet")
    return k, len(D), len(R), time.time() - t0


def main():
    times_csv, data, out = sys.argv[1:4]
    workers, hmax, bsize = CP.ephem_workers, E_.ephem_hmax, CP.ephem_batch
    Path(out).mkdir(parents=True, exist_ok=True)
    P = pd.read_csv(times_csv)
    ut = pd.to_datetime(P.ut)
    jd_utc = ut.map(pd.Timestamp.to_julian_date).values
    d0, d1, d2 = E_.delta_t
    dT = (d0 + d1 * (ut.dt.year + ut.dt.day_of_year / 365.25 - d2)).values / 86400      # Delta T (TT - UT), linear near 1950
    rows = []
    for when, off in [("start", -0.5), ("mid", 0.0), ("end", 0.5)]:
        j = jd_utc + off * P.exp.values / 1440
        rows.append(pd.DataFrame({"plate": P.plate, "when": when, "jd_utc": j, "tdb": j + dT, "ra_c": P.ra_c, "dec_c": P.dec_c}))
    times = pd.concat(rows, ignore_index=True)
    t0 = time.time(); D = read_mpcorb(Path(data) / "MPCORB.DAT.gz", hmax)
    if os.environ.get("ONLY_DESIG"):                                   # only these orbits (those already known near the fields)
        keep = set(Path(os.environ["ONLY_DESIG"]).read_text().split()); D = D[D.desig.isin(keep)]
    print(f"{len(D):,} orbits with H <= {hmax} read in {time.time() - t0:.0f} s; {len(times)} times", flush=True)
    D = D.sort_values("a").reset_index(drop=True)                      # similar orbits share a batch (similar step sizes)
    batches = [(k, D.iloc[s:s + bsize], times, data, out) for k, s in enumerate(range(0, len(D), bsize))]
    done = {int(f.stem[-4:]) for f in Path(out).glob("asteroids_part*.parquet")}
    with ProcessPoolExecutor(workers) as ex:
        for k, n, m, sec in ex.map(run_batch, [b for b in batches if b[0] not in done]):
            print(f"batch {k}: {n} orbits, {m} in-field positions, {sec:.0f} s", flush=True)


if __name__ == "__main__":
    main()
