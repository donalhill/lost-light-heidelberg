"""The pass decision, the cut-flow, and the open record of the selection (notes/requirements.md D2, D7, B2). Pure
functions of a flags table (one boolean column fail_<cut> per cut in selection.CUTS) so that every rule is tested once."""

from heidelberg.lib.selection import CUTS, DECISION

# The fields of a pair's open selection record in blind mode (B2): nothing that counts or describes the candidates.
OPEN_FIELDS = ("pair", "epoch", "chance_final", "chance_dirs", "injected", "chance_by_radius", "injections_passed",
               "G_eff_half", "cutflow", "source", "vetoes_complete", "offline", "rule", "fidelity_source",
               "asteroid_ephemeris", "scanner_marks", "cache_used", "run_key", "blind", "refused")


def passed(T):
    """True where a row fails none of the deciding cuts."""
    return ~T[[f"fail_{c}" for c in DECISION]].any(axis=1)


def cutflow(T):
    """Per population: rows, passed, remote_evaluated, and per cut [fails it, fails it alone]: alone means it is the
    only cut the row fails."""
    flow = {}
    for pop, g in T.groupby("population"):
        nd = g[[f"fail_{c}" for c in DECISION]].sum(axis=1).to_numpy()
        d = {}
        for c in CUTS:
            f = g[f"fail_{c}"].to_numpy(bool)
            d[c] = [int(f.sum()), int((f & (nd == 1)).sum())]
        d.update(rows=int(len(g)), passed=int(passed(g).sum()), remote_evaluated=int(g["remote_evaluated"].sum()) if "remote_evaluated" in g else 0)
        flow[pop] = d
    return flow


def open_record(r):
    """A pair's selection record as written to the open outputs in blind mode: the whitelisted fields, with the
    candidates' cut-flow removed."""
    out = {k: v for k, v in r.items() if k in OPEN_FIELDS}
    if "cutflow" in out:
        out["cutflow"] = {k: v for k, v in out["cutflow"].items() if k != "candidate"}
    out["blind"] = True
    return out
