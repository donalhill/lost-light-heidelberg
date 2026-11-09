"""The background-constrained likelihood for the shadow test (notes/pipeline_spec.md). Per pair i a background B_i,
constrained by the pair's chance estimate c_i (a mean over ns shifts) through the quasi-Poisson term (ns / OD)(c_i ln B_i
- B_i), and glints of density eta suppressed in shadow by kappa. For finals with shadow shares s_ij:
    -lnL = sum_i [ B_i + eta E_i (1 - kappa sbar_i) - sum_j ln(B_i w_ij + eta E_i (1 - kappa s_ij)) - (ns / OD)(c_i ln B_i - B_i) ],
with sbar_i the mean shadow share of the pair's controls, w_ij the background's map at each final (mean 1 over the
pair's controls; default 1), and E_i the pair's efficiency x exposure (default 1). The test is kappa > 0 with eta and
every B_i profiled: q = 2 (lnL(kappa hat) - lnL(0)) for kappa hat > 0, else 0, its null distribution from toys.
The background's level is free (#53): every fit first multiplies the chance estimates by one common factor, the number
of surviving sources over their total (scale_hat), so that they set only how the background varies from pair to pair
and the data set its level. A constant factor on every chance estimate thus changes no fit, statistic, or limit; the
toys repeat the step, so the calibration includes it."""
import numpy as np
from numba import njit
from scipy.optimize import minimize, minimize_scalar

from heidelberg.config import A, DERIVED

T_ = A.statistics
NS = A.background.shifts                                   # shift directions in the chance estimate
B_LO, B_HI = T_.b_lo, T_.b_hi                              # bounds on each pair's background
OD = DERIVED["overdispersion"]                      # the chance estimate's variance over its Poisson variance (measured on the tuning pairs)
EXCESS = OD - 1                                            # the variance's excess over Poisson, for the toys


def pad(fin):
    """Finals' shares per pair (list of arrays) as a padded matrix and mask."""
    n = max(1, max((len(x) for x in fin), default=1)); S = np.zeros((len(fin), n)); m = np.zeros((len(fin), n), bool)
    for i, x in enumerate(fin):
        S[i, :len(x)] = x; m[i, :len(x)] = True
    return S, m


def fit_gradient(P, C):
    """The background's map in box coordinates, pooled over pairs: density proportional to exp(g . r) over each pair's
    valid area (represented by its controls). P, C: per pair arrays (n, 2) of the chance passes' and the controls'
    coordinates (east and north offsets over the box's half-widths). Returns g (2,) by maximum likelihood."""
    def nll(g):
        return -sum((p @ g).sum() - len(p) * np.log(np.mean(np.exp(c @ g))) for p, c in zip(P, C) if len(p) and len(c))
    return minimize(nll, [0.0, 0.0], method="Nelder-Mead", options=dict(xatol=T_.grad_xatol, fatol=T_.grad_fatol)).x


@njit(cache=True)
def _block(a, lo, n):
    """numpy's pairwise_sum of a[lo:lo + n] for n <= 128: in sequence below 8 terms, else 8 running sums combined in pairs
    and the remainder added in sequence."""
    if n < 8:
        r = 0.0
        for i in range(lo, lo + n):
            r += a[i]
        return r
    r0, r1, r2, r3, r4, r5, r6, r7 = a[lo], a[lo + 1], a[lo + 2], a[lo + 3], a[lo + 4], a[lo + 5], a[lo + 6], a[lo + 7]
    i = 8
    while i < n - n % 8:
        j = lo + i
        r0 += a[j]; r1 += a[j + 1]; r2 += a[j + 2]; r3 += a[j + 3]; r4 += a[j + 4]; r5 += a[j + 5]; r6 += a[j + 6]; r7 += a[j + 7]
        i += 8
    res = ((r0 + r1) + (r2 + r3)) + ((r4 + r5) + (r6 + r7))
    while i < n:
        res += a[lo + i]; i += 1
    return res


@njit(cache=True)
def _pairwise(a, lo, n):
    """The sum of a[lo:lo + n] in numpy's order (pairwise_sum: blocks of at most 128 by _block, longer runs split in halves
    at a multiple of 8 and the halves added), so that the compiled step adds a pair's terms exactly as numpy does. The
    halving is walked with an explicit stack (numba's cache does not load a recursive function safely)."""
    if n <= 128:
        return _block(a, lo, n)
    L = np.empty(64, np.int64); N = np.empty(64, np.int64); ST = np.zeros(64, np.int64); LEFT = np.empty(64)
    sp = 0; L[0] = lo; N[0] = n; ret = 0.0
    while sp >= 0:
        k = N[sp]
        if k <= 128:                              # a block: its sum goes to the node above
            ret = _block(a, L[sp], k); sp -= 1
            continue
        n2 = k // 2; n2 -= n2 % 8
        if ST[sp] == 0:                           # first the left half
            ST[sp] = 1; sp += 1; L[sp] = L[sp - 1]; N[sp] = n2; ST[sp] = 0
        elif ST[sp] == 1:                         # the left half's sum is in ret: then the right half
            LEFT[sp] = ret; ST[sp] = 2; sp += 1; L[sp] = L[sp - 1] + n2; N[sp] = k - n2; ST[sp] = 0
        else:                                     # both halves done
            ret = LEFT[sp] + ret; sp -= 1
    return ret


@njit(cache=True)
def _newton_step(B, x, c, Wm, g, m, k, lo, hi, maxstep):
    """One Newton step in ln B_i for every pair, given B = exp(x): nll_profiled_ref's arithmetic in numpy's order (each
    pair's padded row of terms, masked entries zero, summed by _pairwise). No exp or log here: numpy computes those, so
    the two forms are the same operations on any machine. Returns the new x and the largest move."""
    M, n = Wm.shape
    xn = np.empty(M); t1 = np.empty(n); t2 = np.empty(n); dmax = 0.0
    for i in range(M):
        for j in range(n):
            inv = 1.0 / max(B[i] * Wm[i, j] + g[i, j], 1e-300) if m[i, j] else 0.0
            t1[j] = Wm[i, j] * inv; t2[j] = Wm[i, j] * Wm[i, j] * inv * inv
        S1 = _pairwise(t1, 0, n); S2 = _pairwise(t2, 0, n)
        f1 = B[i] * S1 - B[i] + k * (c[i] - B[i]); f2 = B[i] * S1 - B[i] - k * B[i] - B[i] * B[i] * S2
        step = -f1 / f2 if f2 < 0 else np.sign(f1)
        xn[i] = min(max(x[i] + min(max(step, -maxstep), maxstep), lo), hi); dmax = max(dmax, abs(xn[i] - x[i]))
    return xn, dmax


def nll_profiled(eta, kappa, c, sbar, S, m, od, E=None, ns=NS, W=None):
    """-lnL at (eta, kappa) with every B_i at its profile maximum (Newton's method in ln B_i), its inner step compiled;
    identical to nll_profiled_ref (T9)."""
    if eta < 0 or np.any((1 - kappa * S)[m] < 0) or np.any(1 - kappa * sbar < 0):
        return 1e12
    E = np.ones(len(sbar)) if E is None else np.asarray(E, float)
    g = np.where(m, eta * E[:, None] * (1 - kappa * S), 0.0); Wm = np.where(m, 1.0 if W is None else W, 0.0)
    k = ns / od; c = np.asarray(c, float); lo, hi = np.log(B_LO), np.log(B_HI)
    x = np.clip(np.log(np.maximum(c, T_.fit_eta_start)), lo, hi)
    for _ in range(T_.newton_iter):
        xn, dmax = _newton_step(np.exp(x), x, c, Wm, g, m, k, lo, hi, T_.newton_step)
        if dmax < T_.newton_tol:
            x = xn; break
        x = xn
    B = np.exp(x); d = B[:, None] * Wm + g
    ll = np.where(m, np.log(np.maximum(d, 1e-300)), 0.0).sum(1) - B - eta * E * (1 - kappa * sbar) + k * (c * np.log(B) - B)
    return -float(ll.sum())


def nll_profiled_ref(eta, kappa, c, sbar, S, m, od, E=None, ns=NS, W=None):
    """-lnL at (eta, kappa) with every B_i at its profile maximum, by Newton's method in ln B_i, vectorised over pairs: the
    reference form, in numpy, against which the compiled nll_profiled is tested."""
    if eta < 0 or np.any((1 - kappa * S)[m] < 0) or np.any(1 - kappa * sbar < 0):
        return 1e12
    E = np.ones(len(sbar)) if E is None else np.asarray(E, float)
    g = np.where(m, eta * E[:, None] * (1 - kappa * S), 0.0); Wm = np.where(m, 1.0 if W is None else W, 0.0)
    k = ns / od; c = np.asarray(c, float); lo, hi = np.log(B_LO), np.log(B_HI)
    x = np.clip(np.log(np.maximum(c, T_.fit_eta_start)), lo, hi)
    for _ in range(T_.newton_iter):
        B = np.exp(x); d = B[:, None] * Wm + g; inv = np.where(m, 1.0 / np.maximum(d, 1e-300), 0.0)
        S1 = (Wm * inv).sum(1); S2 = (Wm * Wm * inv * inv).sum(1)
        f1 = B * S1 - B + k * (c - B)                          # d lnL_i / d ln B
        f2 = B * S1 - B - k * B - B * B * S2                   # d2 lnL_i / d (ln B)2
        step = np.where(f2 < 0, -f1 / np.where(f2 < 0, f2, -1.0), np.sign(f1))
        xn = np.clip(x + np.clip(step, -T_.newton_step, T_.newton_step), lo, hi)
        if np.max(np.abs(xn - x)) < T_.newton_tol:
            x = xn; break
        x = xn
    B = np.exp(x); d = B[:, None] * Wm + g
    ll = np.where(m, np.log(np.maximum(d, 1e-300)), 0.0).sum(1) - B - eta * E * (1 - kappa * sbar) + k * (c * np.log(B) - B)
    return -float(ll.sum())


def scale_hat(c, fin):
    """The chance estimates' common factor: the number of surviving sources over the chance estimates' total (1 if either
    is zero), reported with each fit; near 1 when the chance estimates' total agrees with the surviving sources'."""
    n = sum(len(x) for x in fin); C = float(np.sum(c))
    return n / C if n > 0 and C > 0 else 1.0


def fit(c, sbar, fin, od=OD, kappa_fixed=None, E=None, ns=NS, w=None):
    """(minimum -lnL, (ln eta, kappa)) with kappa free in [kappa_lo, kappa_hi] or fixed, on the rescaled chance estimates."""
    c = np.asarray(c, float) * scale_hat(c, fin)
    S, m = pad(fin); W = None if w is None else pad([np.asarray(x, float) for x in w])[0] + ~m
    if kappa_fixed is not None:
        r = minimize(lambda x: nll_profiled(np.exp(x[0]), kappa_fixed, c, sbar, S, m, od, E, ns, W), [np.log(T_.fit_eta_start)], method="Nelder-Mead",
                     options=dict(xatol=T_.fit_xatol, fatol=T_.fit_fatol))
        return r.fun, np.r_[r.x, kappa_fixed]
    best = None
    for k0 in T_.fit_kappa_starts:
        r = minimize(lambda x: nll_profiled(np.exp(x[0]), x[1], c, sbar, S, m, od, E, ns, W) if T_.kappa_lo <= x[1] <= T_.kappa_hi else 1e12, [np.log(T_.fit_eta_start), k0],
                     method="Nelder-Mead", options=dict(xatol=T_.fit_xatol, fatol=T_.fit_fatol))
        if best is None or r.fun < best.fun:
            best = r
    return best.fun, best.x


def q_stat(c, sbar, fin, od=OD, E=None, ns=NS, w=None):
    """The test statistic and kappa hat (w: per pair, the background weight at each final)."""
    f0, _ = fit(c, sbar, fin, od, 0.0, E, ns, w); f1, x1 = fit(c, sbar, fin, od, None, E, ns, w)
    return (2 * (f0 - f1) if x1[1] > 0 else 0.0), float(x1[1])


def ul_qt(c, sbar, fin, f0, od=OD, E=None, ns=NS, w=None, Bref=None, scale=False):
    """The limit's statistics at each glint share in f0 (kappa = 1, every background profiled), from the same fits: the
    one-sided q = 2 (-lnL(eta0) + lnL(eta hat)), 0 where eta hat lies above eta0, and the two-sided t, the same ratio
    without the zero (eta hat not negative: the ordering of Feldman and Cousins), on the rescaled chance estimates;
    eta0 = Bref f0 / (1 - f0), Bref a fixed reference background (default mean(c); the stage passes the sample's mean
    chance estimate and the calibration the tuning pairs' background, its forecast of that; a constant, not rescaled by
    the fit, so that a share means one glint density whatever the chance estimates' scale). With q_crit(f0) from toys at each f0 (a Neyman belt), the limit is the largest f0 with q <= q_crit; with
    t_crit(f0), the interval is every f0 with t <= t_crit. Returns an array of shape (2, len(f0)), and with scale also
    the factor on the chance estimates (scale_hat)."""
    s = scale_hat(c, fin); Bref = float(np.mean(c)) if Bref is None else float(Bref); c = np.asarray(c, float) * s
    S, m = pad(fin); W = None if w is None else pad([np.asarray(x, float) for x in w])[0] + ~m
    nll = lambda e: nll_profiled(e, 1.0, c, sbar, S, m, od, E, ns, W)
    e0 = T_.ul_eta_min; Bm = max(Bref, e0); emax = T_.ul_eta_max_factor * Bm
    r = minimize_scalar(nll, bounds=(e0, emax), method="bounded", options=dict(xatol=T_.ul_xatol))
    e_hat, n_min = (r.x, r.fun) if r.fun < nll(e0) else (0.0, nll(e0))
    t = [(e, max(2 * (nll(e) - n_min), 0.0)) for e in Bm * np.asarray(f0) / (1 - np.asarray(f0))]
    qt = np.array([[0.0 if e_hat >= e else x for e, x in t], [x for _, x in t]])
    return (qt, s) if scale else qt


def ul_q(c, sbar, fin, f0, od=OD, E=None, ns=NS, w=None, Bref=None):
    """The one-sided statistic alone (ul_qt)."""
    return ul_qt(c, sbar, fin, f0, od, E, ns, w, Bref)[0]
