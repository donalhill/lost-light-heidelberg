"""Cross-matches against the primary archives (CDS is not used: R8): Gaia DR3 from
GAVO and IRSA (TAP joins on an uploaded table of positions), Pan-STARRS DR2 from MAST, and VSX, the cataclysmic binaries, and the
Asiago supernovae from the frozen local snapshots that heidelberg/stages/snapshots.py writes (SNAP, set by
selection.configure). Returns what CDS XMatch returns for these catalogues: a table of all pairs within the radius,
with column i (the position's row) and the VizieR column names selection.py reads. A catalogue not served here, or a
position outside every VSX field downloaded, is an error (never a silent "no match")."""
import io
import time

import numpy as np
import requests
from astropy.table import Table, vstack

from heidelberg.config import A, C

N_, W_ = C.network, C.processes

IRSA = "https://irsa.ipac.caltech.edu/TAP"
GAVO = "https://dc.g-vo.org/tap"
ESA = "https://gea.esac.esa.int/tap-server/tap"
MAST_PS1 = "https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/dr2/mean.csv"
GAIA_TAP = ((ESA, "gaiadr3.gaia_source"), (IRSA, "gaia_dr3_source"))   # one DR3 table: ESA's, else IRSA's copy (identical rows)
DEAD = set()                                                         # services that failed in this process: not asked again


def host_snapshots(kmax):
    """The supernova hosts' frozen all-sky lists (heidelberg/stages/snapshots.py writes them; selection.sn_hosts reads
    them): file name -> (IRSA query, fewest rows of a complete answer), for 2MASS extended sources and the Large Galaxy
    Atlas with K < kmax. One query each with no position constraint (cones covering the whole field set missed about
    1.5% of the sky's galaxies)."""
    tag = f"k{round(kmax * 10)}"
    return {f"tmass_xsc_{tag}.ecsv": (f"SELECT ra, dec, r_k20fe FROM fp_xsc WHERE k_m_k20fe < {kmax}", 131000),
            f"lga_v2_{tag}.ecsv": (f"SELECT ra, dec, r_k20fe FROM lga_v2 WHERE k_m_k20fe < {kmax}", 600)}


class Timed(requests.Session):
    """A session whose requests fail instead of hanging: connect_timeout_s to connect, read_timeout_s to answer
    (config/compute.json; requirement R8)."""
    def request(self, *a, **kw):
        kw.setdefault("timeout", (N_.connect_timeout_s, N_.read_timeout_s)); return super().request(*a, **kw)


def gaia_async(cols, where, maxrec=N_.async_maxrec, tries=N_.async_tries, wait=N_.async_wait_s):
    """Rows of Gaia DR3 gaia_source (cols, a WHERE clause) by an async TAP job at ESA, else at IRSA's copy of the same
    table, each with a timed session and retried; a service that fails is not asked again in this process."""
    import pyvo
    last = None
    for url, table in GAIA_TAP:
        if url in DEAD:
            continue
        svc = pyvo.dal.TAPService(url, session=Timed())
        for k in range(tries):
            try:
                j = svc.submit_job(f"SELECT {cols} FROM {table} WHERE {where}", maxrec=maxrec); j.run()
                j.wait(phases=["COMPLETED", "ERROR", "ABORTED"], timeout=wait)
                if j.phase != "COMPLETED":
                    raise RuntimeError(f"{url}: {j.phase}")
                t = j.fetch_result().to_table(); j.delete(); return t
            except Exception as e:
                last = e
                if k < tries - 1:
                    time.sleep(N_.tap_backoff_s * (k + 1))
        DEAD.add(url)
    raise RuntimeError(f"no Gaia service answered: {last}")


def _tap_upload(url, adql, ra, dec, chunk=N_.tap_chunk, tries=N_.tap_tries):
    """Synchronous TAP queries joining an uploaded table of positions (idx, ra_in, dec_in), in chunks, retried."""
    import pyvo
    svc = pyvo.dal.TAPService(url, session=Timed()); out = []
    for a in range(0, len(ra), chunk):
        up = Table({"idx": np.arange(a, min(a + chunk, len(ra))), "ra_in": np.asarray(ra, float)[a:a + chunk], "dec_in": np.asarray(dec, float)[a:a + chunk]})
        for k in range(tries):
            try:
                r_ = svc.run_sync(adql, uploads={"t": up}, maxrec=N_.tap_maxrec).to_table(); r_.rename_column("idx", "i"); out.append(r_); break
            except Exception:
                if k == tries - 1:
                    raise
                time.sleep(N_.tap_backoff_s * (k + 1))
    out = [t for t in out if len(t)]
    return vstack(out) if out else None


def gaia(ra, dec, radius):
    """Gaia DR3 by an upload join: GAVO's gaia.dr3lite first, then IRSA's gaia_dr3_source (the ESA archive scans the
    whole table for upload joins and aborts). The two give identical answers (at 21,964 positions, 19,094
    pairs, G equal to 1e-6 mag)."""
    qg = (f"SELECT u.idx, g.ra, g.dec, g.pmra, g.pmdec, g.phot_g_mean_mag FROM TAP_UPLOAD.t AS u JOIN gaia.dr3lite AS g "
          f"ON DISTANCE(POINT(g.ra, g.dec), POINT(u.ra_in, u.dec_in)) < {radius / 3600.0}")
    q = (f"SELECT u.idx, g.ra, g.dec, g.pmra, g.pmdec, g.phot_g_mean_mag FROM TAP_UPLOAD.t u, gaia_dr3_source g "
         f"WHERE CONTAINS(POINT('ICRS', g.ra, g.dec), CIRCLE('ICRS', u.ra_in, u.dec_in, {radius / 3600.0})) = 1")
    try:
        m = _tap_upload(GAVO, qg, ra, dec, tries=N_.gavo_tries); svc = GAVO
    except Exception:
        m = _tap_upload(IRSA, q, ra, dec); svc = IRSA
    if m is not None:
        m.meta["service"] = svc
        for a, b in (("ra", "RAdeg"), ("dec", "DEdeg"), ("pmra", "pmRA"), ("pmdec", "pmDE"), ("phot_g_mean_mag", "Gmag")):
            m.rename_column(a, b)
    return m


def _cones(fn, ra, dec, radius, workers=W_.cone_workers):
    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(workers) as ex:
        res = list(ex.map(lambda k: fn(k, ra[k], dec[k], radius), range(len(ra))))
    res = [t for t in res if t is not None and len(t)]
    return vstack(res) if res else None


def _get(url, params, tries=N_.http_tries):
    for k in range(tries):
        try:
            r = requests.get(url, params=params, timeout=N_.http_timeout_s)
            if r.ok:
                return r.text
        except requests.exceptions.RequestException:
            pass
        time.sleep(N_.http_backoff_s * (k + 1))
    raise RuntimeError(f"{url} failed")


def ps1(ra, dec, radius):
    def one(k, r, d, rad):
        txt = _get(MAST_PS1, {"ra": r, "dec": d, "radius": rad / 3600.0, "nDetections.gte": A.selection.ps1_min_detections, "columns": "[raMean,decMean,gMeanPSFMag,gMeanKronMag]"})
        import pandas as pd
        t = pd.read_csv(io.StringIO(txt)) if txt.strip().count("\n") >= 1 else None
        if t is None or not len(t):
            return None
        g, gk = (np.where(t[c].to_numpy(float) <= -999, np.nan, t[c].to_numpy(float)) for c in ("gMeanPSFMag", "gMeanKronMag"))
        sep = np.degrees(np.arccos(np.clip(_unit(t.raMean, t.decMean) @ _unit([r], [d])[0], -1, 1))) * 3600   # as CDS XMatch's angDist
        return Table({"i": np.full(len(t), k), "gmag": g, "gKmag": gk, "angDist": sep})
    m = _cones(one, np.asarray(ra, float), np.asarray(dec, float), radius)
    if m is not None:
        m.meta["service"] = MAST_PS1
    return m


PS1_HATS = "stpubdata/panstarrs/ps1/public/hats/otmo/otmo"   # Pan-STARRS DR2 ObjectThin + MeanObject, HATS (release 2025-08-22), AWS open data
PS1_COLS = ["raMean", "decMean", "gMeanPSFMag", "gMeanKronMag", "nDetections"]


def ps1_hats(ra, dec, radius, workers=W_.ps1_workers):
    """Pan-STARRS DR2 mean objects from STScI's HATS copy on S3 (the table MAST's mean API serves), as ps1() returns
    them: every object with nDetections >= 1 within radius of each position (i, gmag, gKmag, angDist). The tiles read are
    those any position's disc touches (healpy query_disc, inclusive), so positions near a tile edge see the neighbours."""
    import healpy as hp
    import pandas as pd
    import pyarrow.fs as pf
    import pyarrow.parquet as pq
    from concurrent.futures import ThreadPoolExecutor
    from scipy.spatial import cKDTree
    ra, dec = np.asarray(ra, float), np.asarray(dec, float); U = _unit(ra, dec); rr = np.radians(radius / 3600)
    fs = pf.S3FileSystem(anonymous=True, region="us-east-1")
    if "parts" not in SNAP:
        with fs.open_input_file(f"{PS1_HATS}/partition_info.csv") as fh:
            SNAP["parts"] = set(map(tuple, pd.read_csv(fh)[["Norder", "Npix"]].to_numpy()))
    parts = SNAP["parts"]; need = set()
    for o in sorted({o for o, _ in parts}):
        ns = hp.order2nside(o)
        for u in U:
            need |= {(o, int(p)) for p in hp.query_disc(ns, u, rr, inclusive=True, nest=True) if (o, int(p)) in parts}
    tu = cKDTree(U); chord = 2 * np.sin(rr / 2)
    def read(op):                                   # each tile cut at once to the objects near some position
        o, p = op
        t = pq.read_table(f"{PS1_HATS}/dataset/Norder={o}/Dir={p // 10000 * 10000}/Npix={p}.parquet", filesystem=fs, columns=PS1_COLS).to_pandas()
        d_, _ = tu.query(_unit(t.raMean.to_numpy(float), t.decMean.to_numpy(float)), distance_upper_bound=chord * 1.001)
        return t[np.isfinite(d_)]
    with ThreadPoolExecutor(workers) as ex:
        T = pd.concat(list(ex.map(read, sorted(need))), ignore_index=True) if need else pd.DataFrame(columns=PS1_COLS)
    T = T[T.nDetections >= A.selection.ps1_min_detections]
    if not len(T):
        return None
    V = _unit(T.raMean.to_numpy(float), T.decMean.to_numpy(float))
    lst = cKDTree(V).query_ball_point(U, 2 * np.sin(rr / 2))
    i = np.array([k for k, js in enumerate(lst) for _ in js], int); j = np.array([jj for js in lst for jj in js], int)
    if not len(i):
        return None
    g, gk = (np.where(T[c].to_numpy(float)[j] <= -999, np.nan, T[c].to_numpy(float)[j]) for c in ("gMeanPSFMag", "gMeanKronMag"))
    sep = np.degrees(np.arccos(np.clip(np.einsum("ij,ij->i", U[i], V[j]), -1, 1))) * 3600
    out = Table({"i": i, "gmag": g, "gKmag": gk, "angDist": sep}); out.meta["service"] = f"s3://{PS1_HATS}"; return out


SNAP = {"dir": None}                                      # the snapshot directory (selection.configure(snapshots=...))
LOCAL = ("B/vsx/vsx", "B/cb/cbdata", "B/sn/sncat")


def _unit(r, d):
    r, d = np.radians(np.asarray(r, float)), np.radians(np.asarray(d, float)); return np.c_[np.cos(d) * np.cos(r), np.cos(d) * np.sin(r), np.sin(d)]


def _within(ra, dec, sra, sde, radius):
    """All (position, source) pairs within radius arcsec, as a table with column i (the position's row)."""
    from scipy.spatial import cKDTree
    if not len(sra):
        return None
    lst = cKDTree(_unit(sra, sde)).query_ball_point(_unit(ra, dec), 2 * np.sin(np.radians(radius / 3600) / 2))
    i = np.array([k for k, js in enumerate(lst) for _ in js], int)
    return Table({"i": i}) if len(i) else None


def snapshot_files():
    import glob
    import os
    d = SNAP["dir"]
    if not d or not os.path.isdir(d):
        raise RuntimeError("catalogue snapshots not set or missing (heidelberg/stages/snapshots.py; HD_SNAPSHOTS)")
    vsx = [f"{d}/B_vsx_vsx.tsv"] if os.path.exists(f"{d}/B_vsx_vsx.tsv") else sorted(glob.glob(f"{d}/vsx/*.json"))   # the whole catalogue, else the fields
    return vsx + [f"{d}/B_cb_cbdata.tsv", f"{d}/B_sn_sncat.tsv"]


def vsx_local(ra, dec, radius):
    """VSX within radius arcsec: from the whole catalogue (SNAP_DIR/B_vsx_vsx.tsv, VizieR B/vsx/vsx, fetched once by the
    Snakefile's fetch_vsx) when present, else from the per-field cone answers (SNAP_DIR/vsx/*.json), every position
    inside a downloaded field."""
    import json
    import os
    ra, dec = np.asarray(ra, float), np.asarray(dec, float)
    if os.path.exists(f"{SNAP['dir']}/B_vsx_vsx.tsv"):
        return _within(ra, dec, *_vizier_tsv("B/vsx/vsx"), radius)
    fs = [f for f in snapshot_files() if f.endswith(".json")]
    cover = np.zeros(len(ra), bool); seen = {}
    for f in fs:
        js = json.load(open(f)); fd = js["field"]
        cover |= np.degrees(np.arccos(np.clip(_unit(ra, dec) @ _unit([fd["ra"]], [fd["dec"]])[0], -1, 1))) < fd["radius_deg"] - radius / 3600
        for o in js["objects"]:
            seen[o.get("Name")] = (float(o["RA2000"]), float(o["Declination2000"]))
    if not cover.all():
        raise RuntimeError(f"VSX: {int((~cover).sum())} of {len(ra)} positions outside every downloaded field")
    v = np.array(list(seen.values()), float).reshape(-1, 2)
    return _within(ra, dec, v[:, 0], v[:, 1], radius)


def _vizier_tsv(cat):
    import pandas as pd
    f = f"{SNAP['dir']}/{cat.replace('/', '_')}.tsv"; snapshot_files()
    T = pd.read_csv(f, sep="\t", comment="#", skip_blank_lines=True, dtype=str).iloc[2:]     # after the units and dashes rows
    r, d = pd.to_numeric(T["_RAJ2000"], errors="coerce").to_numpy(), pd.to_numeric(T["_DEJ2000"], errors="coerce").to_numpy()
    ok = np.isfinite(r) & np.isfinite(d); return r[ok], d[ok]


def cb_local(ra, dec, radius):
    return _within(np.asarray(ra, float), np.asarray(dec, float), *_vizier_tsv("B/cb/cbdata"), radius)


def sn_local(ra, dec, radius):
    return _within(np.asarray(ra, float), np.asarray(dec, float), *_vizier_tsv("B/sn/sncat"), radius)


PS1_SOURCE = "hats"                                       # set by heidelberg/stages/select.py (--ps1)


def ps1_any(ra, dec, radius):
    """Pan-STARRS from the HATS copy on S3 (PS1_SOURCE hats, fast inside AWS), else MAST's mean API (one call per
    position). Identical answers (at 150 positions, 332 pairs; magnitudes to the 0.001 of MAST's CSV)."""
    return ps1_hats(ra, dec, radius) if PS1_SOURCE == "hats" else ps1(ra, dec, radius)


SERVED = {"I/355/gaiadr3": gaia, "II/349/ps1": ps1_any, "B/vsx/vsx": vsx_local, "B/cb/cbdata": cb_local, "B/sn/sncat": sn_local}


def xmatch(ra, dec, cat, radius):
    if not len(ra):
        return None
    if cat not in SERVED:
        raise RuntimeError(f"{cat}: not served by the primary archives or the snapshots")
    return SERVED[cat](ra, dec, radius)
