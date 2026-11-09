"""A pair's efficiency E_i (notes/requirements.md E1): the share of star-like injections (duration 0, no sharpening)
introduced on both plates that pass every deciding cut. One definition, used by the shadow stage and the role report."""


def pair_efficiency(J, flags, g_max=None):
    """J: the pair's injection table (PAIR_inj.parquet); flags: its flags table (PAIR_flags.parquet), whose injection
    rows carry the injection's index in column inj and the verdict in column passed. g_max: only flashes brighter than
    this G (the rate's brighter range)."""
    star = J[(J.dur == 0) & (J.sharpen == 1) & J.placed_a.astype(bool) & J.placed_b.astype(bool)] if "dur" in J else J
    if g_max is not None:
        star = star[star.G < g_max]
    inj = flags[flags.population == "injection"]
    ok = inj[inj.passed].inj if "inj" in inj else []
    return float(star.index.isin(ok).mean()) if len(star) else 0.0
