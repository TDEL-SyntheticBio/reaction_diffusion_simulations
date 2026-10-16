"""
One estimator for point-pattern order, used identically on in vitro wells and
simulation fields.

g(r):  border-method pair correlation. Raw pair counts and the CSR denominator
       are returned separately so several windows can be pooled exactly:
           g_pooled = sum(counts) / sum(denom)
       In reduced units x = r/a with a = 1/sqrt(lambda), the denominator is
       density-free, so wells of different density pool without smearing.

S(k):  exact finite-window structure factor of the centroids on the reciprocal
       lattice of the window, S(k) = |sum_j exp(-i k . r_j)|^2 / N.
       CSR reference is S = 1 at every non-zero k.

Nulls: CSR            - uniform points, no interaction, no hard core
       hard-core RSA  - non-overlapping discs at the observed diameters
       germ-grain     - discs at the observed diameters placed independently
                        WITH overlap, then re-segmented as connected components.
                        This is the null that reproduces "random nucleation plus
                        merging" and is the one that matters for the hole.
"""
import numpy as np
from scipy.spatial import cKDTree
from scipy import ndimage as ndi


# ----------------------------------------------------------------- g(r)

def gr_counts(cy, cx, H, W, edges):
    """Border-method pair counts in the annuli given by `edges` (same length
    units as cy/cx/H/W). Returns (counts, denom, n_focal, lam) or None."""
    rmax = float(edges[-1])
    pts = np.column_stack([np.asarray(cy, float), np.asarray(cx, float)])
    n = len(pts)
    if n < 10:
        return None
    lam = n / (H * W)
    inner = ((pts[:, 0] >= rmax) & (pts[:, 0] <= H - rmax) &
             (pts[:, 1] >= rmax) & (pts[:, 1] <= W - rmax))
    foc = pts[inner]
    if len(foc) < 10:
        return None
    tree = cKDTree(pts)
    d = tree.sparse_distance_matrix(cKDTree(foc), rmax, output_type='ndarray')
    dd = d['v']
    dd = dd[dd > 1e-9]
    counts = np.histogram(dd, bins=edges)[0].astype(float)
    ring = np.pi * (np.asarray(edges[1:]) ** 2 - np.asarray(edges[:-1]) ** 2)
    denom = len(foc) * lam * ring
    return counts, denom, len(foc), lam


def gr_reduced(cy, cx, H, W, x_edges):
    """g(r) on a grid of reduced separations x = r/a, a = 1/sqrt(lambda).
    The returned denominator is density-free, so these pool across wells."""
    n = len(cy)
    if n < 10:
        return None
    lam = n / (H * W)
    a = 1.0 / np.sqrt(lam)
    out = gr_counts(cy, cx, H, W, np.asarray(x_edges) * a)
    if out is None:
        return None
    counts, denom, nfoc, _ = out
    return counts, denom, nfoc, a


def pool(results):
    """Pool a list of (counts, denom, ...) into g with Poisson errors."""
    C = np.sum([r[0] for r in results], axis=0)
    D = np.sum([r[1] for r in results], axis=0)
    g = np.where(D > 0, C / np.maximum(D, 1e-30), np.nan)
    sd = np.where(D > 0, np.sqrt(np.maximum(C, 0)) / np.maximum(D, 1e-30), np.nan)
    return g, sd, C, D


def hole_radius(x_cent, g, C):
    """Largest separation out to which not a single pair was observed."""
    nz = np.nonzero(C > 0)[0]
    if len(nz) == 0:
        return np.nan
    return float(x_cent[nz[0]]) if nz[0] > 0 else 0.0


# ----------------------------------------------------------------- S(k)

def structure_factor(cy, cx, H, W, mmax=90, ka_edges=None):
    """Exact finite-window S(k) binned in reduced wavenumber k*a."""
    y = np.asarray(cy, float); x = np.asarray(cx, float)
    N = len(y)
    if N < 10:
        return None
    lam = N / (H * W)
    a = 1.0 / np.sqrt(lam)
    m = np.arange(-mmax, mmax + 1)
    kx = 2 * np.pi * m / W
    ky = 2 * np.pi * m / H
    A = np.exp(-1j * np.outer(kx, x - x.mean()))
    B = np.exp(-1j * np.outer(ky, y - y.mean()))
    F = B @ A.T
    S = (np.abs(F) ** 2) / N
    KX, KY = np.meshgrid(kx, ky)
    K = np.hypot(KX, KY) * a                      # reduced wavenumber k*a
    ok = K > 1e-12
    if ka_edges is None:
        ka_edges = np.linspace(0, 12, 49)
    idx = np.digitize(K[ok].ravel(), ka_edges) - 1
    nb = len(ka_edges) - 1
    good = (idx >= 0) & (idx < nb)
    cnt = np.bincount(idx[good], minlength=nb).astype(float)
    tot = np.bincount(idx[good], weights=S[ok].ravel()[good], minlength=nb)
    return tot, cnt, a, N


def pool_sk(results):
    """Pool S(k): sum of periodogram values over sum of modes per bin."""
    T = np.sum([r[0] for r in results], axis=0)
    C = np.sum([r[1] for r in results], axis=0)
    S = np.where(C > 0, T / np.maximum(C, 1e-30), np.nan)
    # a flat periodogram has unit-exponential values, so SEM = S/sqrt(modes)
    sem = np.where(C > 0, S / np.sqrt(np.maximum(C, 1)), np.nan)
    return S, sem, C


# ----------------------------------------------------------------- nulls

def csr(n, H, W, rng):
    return rng.uniform(0, H, n), rng.uniform(0, W, n)


def rsa_discs(H, W, diam, rng, tries=4000):
    """Random sequential adsorption of non-overlapping discs, largest first."""
    d = np.sort(np.asarray(diam, float))[::-1]
    cy = np.empty(len(d)); cx = np.empty(len(d)); dd = np.empty(len(d)); p = 0
    for di in d:
        for _ in range(tries):
            yy = rng.uniform(di / 2, H - di / 2)
            xx = rng.uniform(di / 2, W - di / 2)
            if p == 0 or np.all(np.hypot(cy[:p] - yy, cx[:p] - xx) > 0.5 * (dd[:p] + di)):
                cy[p] = yy; cx[p] = xx; dd[p] = di; p += 1
                break
    return cy[:p], cx[:p], dd[:p]


def germ_grain(H, W, diam_pool, phi_target, rng, px=8.0, min_d=40.0, pad=1.0):
    """Boolean model: discs drawn from `diam_pool` placed independently, allowed
    to overlap, rasterised and re-segmented. Germ number is set so the expected
    Boolean coverage equals phi_target: sum(area)/A = -ln(1 - phi).
    Returns (cy, cx, eqd) of the resulting connected components, >= min_d."""
    A = H * W
    need = -np.log(max(1e-9, 1.0 - min(phi_target, 0.95))) * A
    pool_d = np.asarray(diam_pool, float)
    mean_a = np.pi * (pool_d ** 2).mean() / 4.0
    N = max(1, int(round(need / mean_a)))
    d = rng.choice(pool_d, size=N, replace=True)
    # place germ centres in a padded window so border objects are not biased
    cy = rng.uniform(-pad * d.max(), H + pad * d.max(), N)
    cx = rng.uniform(-pad * d.max(), W + pad * d.max(), N)
    ny, nx = int(np.ceil(H / px)), int(np.ceil(W / px))
    m = np.zeros((ny, nx), bool)
    r_px = d / 2.0 / px
    yy_px = cy / px
    xx_px = cx / px
    for k in range(N):
        rr = r_px[k]
        y0 = int(np.floor(yy_px[k] - rr)); y1 = int(np.ceil(yy_px[k] + rr)) + 1
        x0 = int(np.floor(xx_px[k] - rr)); x1 = int(np.ceil(xx_px[k] + rr)) + 1
        y0c, y1c = max(0, y0), min(ny, y1)
        x0c, x1c = max(0, x0), min(nx, x1)
        if y0c >= y1c or x0c >= x1c:
            continue
        gy = np.arange(y0c, y1c)[:, None] + 0.5 - yy_px[k]
        gx = np.arange(x0c, x1c)[None, :] + 0.5 - xx_px[k]
        m[y0c:y1c, x0c:x1c] |= (gy * gy + gx * gx) <= rr * rr
    lab, n = ndi.label(m, structure=np.ones((3, 3), np.uint8))
    if n == 0:
        return np.array([]), np.array([]), np.array([]), float(m.mean())
    ys, xs = np.nonzero(m)
    l = lab[ys, xs]
    cnt = np.bincount(l, minlength=n + 1).astype(float)
    sy = np.bincount(l, weights=ys.astype(float), minlength=n + 1)
    sx = np.bincount(l, weights=xs.astype(float), minlength=n + 1)
    areas = cnt[1:] * px * px
    ccy = (sy[1:] / cnt[1:]) * px
    ccx = (sx[1:] / cnt[1:]) * px
    eqd = 2 * np.sqrt(areas / np.pi)
    keep = eqd >= min_d
    return ccy[keep], ccx[keep], eqd[keep], float(m.mean())


# -------------------------------------------- translation-corrected g(r)

def gr_translation(cy, cx, H, W, edges):
    """Ohser translation-corrected g(r) for a rectangular window. Uses every
    ordered pair, weighting each by A / A_overlap(dx, dy), so no points are
    discarded to a border guard. Returns (num, den, raw_count, n, a)
    with g = sum(num)/sum(den) poolable across windows."""
    pts = np.column_stack([np.asarray(cy, float), np.asarray(cx, float)])
    n = len(pts)
    if n < 10:
        return None
    lam = n / (H * W)
    rmax = float(edges[-1])
    tree = cKDTree(pts)
    d = tree.sparse_distance_matrix(tree, rmax, output_type='ndarray')
    i, j, v = d['i'], d['j'], d['v']
    m = v > 1e-9
    i, j, v = i[m], j[m], v[m]
    if len(v) < 5:
        return None
    dy = np.abs(pts[i, 0] - pts[j, 0])
    dx = np.abs(pts[i, 1] - pts[j, 1])
    Aint = np.maximum((H - dy) * (W - dx), 1e-12)
    w = 1.0 / Aint
    edges = np.asarray(edges, float)
    ring = np.pi * (edges[1:] ** 2 - edges[:-1] ** 2)
    num = np.histogram(v, bins=edges, weights=w)[0] / ring
    cnt = np.histogram(v, bins=edges)[0].astype(float)
    den = np.full(len(ring), lam * lam)
    return num, den, cnt, n, 1.0 / np.sqrt(lam)


def gr_translation_reduced(cy, cx, H, W, x_edges):
    n = len(cy)
    if n < 10:
        return None
    lam = n / (H * W)
    a = 1.0 / np.sqrt(lam)
    out = gr_translation(cy, cx, H, W, np.asarray(x_edges) * a)
    if out is None:
        return None
    num, den, cnt, _, _ = out
    return num, den, cnt, n, a


def pool_translation(results):
    """g with an error bar from the unweighted pair count in each bin."""
    N = np.sum([r[0] for r in results], axis=0)
    D = np.sum([r[1] for r in results], axis=0)
    C = np.sum([r[2] for r in results], axis=0)
    g = np.where(D > 0, N / np.maximum(D, 1e-30), np.nan)
    sd = np.where(C > 0, g / np.sqrt(np.maximum(C, 1)), np.where(D > 0, 0.0, np.nan))
    return g, sd, C
