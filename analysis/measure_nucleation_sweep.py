"""Measure the nucleation sweep with exactly the conventions used everywhere else:
threshold 0.3*b_a, engine even-r hex adjacency, components >= 2 cells (= 41 um
equivalent diameter, matching the in vitro 40 um cut), 28 um cell spacing."""
import os, sys, glob
import numpy as np, pandas as pd
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, 'analysis'))
from analysis import hexgeom
from spatial_order import gr_translation_reduced, pool_translation, structure_factor, pool_sk

UM = 28.0; THR = 1.5; HOUR = 5.4; MIN_CELLS = 2
DAYS = {1: 24/HOUR, 2: 48/HOUR, 3: 72/HOUR, 4: 96/HOUR}
X_E = np.linspace(0, 2.6, 41); X = 0.5*(X_E[1:]+X_E[:-1])
KA = np.linspace(0, 12, 49); KAC = 0.5*(KA[1:]+KA[:-1]); BAND = (KAC >= 1.4) & (KAC <= 3.1)


def measure(a):
    m = a > THR
    ct = hexgeom.component_table(m, min_size=MIN_CELLS)
    ny, nx = a.shape
    area_mm2 = hexgeom.window_area(ny, nx) * UM**2 / 1e6
    n = ct['n_kept']
    if n == 0:
        return dict(n=0, dens=0.0, d_mean=np.nan, d_cv=np.nan, log2sd=np.nan,
                    cov=float(m.mean()), area_mm2=area_mm2)
    areas = ct['area'] * hexgeom.CELL_AREA * UM**2
    eqd = 2*np.sqrt(areas/np.pi)
    return dict(n=n, dens=n/area_mm2, d_mean=float(eqd.mean()),
                d_med=float(np.median(eqd)),
                d_cv=float(eqd.std(ddof=1)/eqd.mean()) if n > 1 else np.nan,
                log2sd=float(np.log2(areas).std(ddof=1)) if n > 1 else np.nan,
                cov=float(m.mean()), area_mm2=area_mm2,
                cy=ct['cy']*UM, cx=ct['cx']*UM, eqd=eqd, ny=ny, nx=nx)


rows, gr_by_rate = [], {}
for f in sorted(glob.glob(os.environ.get('NUCSWEEP_OUT', os.path.join(ROOT, 'results', 'nucleation_sweep')) + '/*.npz')):
    z = np.load(f)
    rate = float(z['nucleation_rate']); rep = int(z['rep']); t = z['t']
    per_day = {}
    for d, tt in DAYS.items():
        k = int(np.argmin(np.abs(t - tt)))
        per_day[d] = measure(z['a'][k].astype(np.float32))
    m4 = per_day[4]
    row = dict(rate=rate, rep=rep, area_mm2=m4['area_mm2'],
               n=m4['n'], dens=m4['dens'], cov=m4['cov'],
               d1=per_day[1]['d_mean'], d2=per_day[2]['d_mean'],
               d3=per_day[3]['d_mean'], d4=m4['d_mean'],
               n1=per_day[1]['n'], n2=per_day[2]['n'], n3=per_day[3]['n'],
               d_cv=m4['d_cv'], log2sd=m4['log2sd'])
    row['g12'] = row['d2']-row['d1']; row['g23'] = row['d3']-row['d2']
    row['g34'] = row['d4']-row['d3']
    if m4['n'] >= 40:
        xmin, xmax, ymin, ymax = hexgeom.window(m4['ny'], m4['nx'])
        H = (ymax-ymin)*UM; W = (xmax-xmin)*UM
        g = gr_translation_reduced(m4['cy'], m4['cx'], H, W, X_E)
        s = structure_factor(m4['cy'], m4['cx'], H, W, ka_edges=KA)
        if g is not None:
            gr_by_rate.setdefault(rate, []).append(g)
            gg, sd, C = pool_translation([g])
            sel = (X >= 0.6) & (X <= 1.4); i1 = int(np.nanargmax(np.where(sel, gg, -9)))
            tail = float(np.nanmean(gg[X >= 1.5]))
            row['amp1'] = gg[i1]-tail
            row['g030'] = gg[int(np.argmin(np.abs(X-0.30)))]
        if s is not None:
            S, _, _ = pool_sk([s]); row['S_band'] = float(np.nanmean(S[BAND]))
    rows.append(row)

df = pd.DataFrame(rows).sort_values(['rate', 'rep'])
df.to_csv(os.path.join(ROOT, 'analysis', 'summaries', 'nucleation_sweep_per_run.csv'), index=False)
agg = df.groupby('rate').agg(['mean', 'std'])
pd.set_option('display.width', 250, 'display.max_columns', 60)

print("="*112)
print("NUCLEATION-RATE SWEEP at the published set (b_a 5, b_i 25, n 10/4, D_i 20, gamma 0.5)")
print("  200x200 = 27.2 mm2, 3 replicates each, read at day 4 (t = 17.8)")
print("="*112)
g = df.groupby('rate')
out = pd.DataFrame({
    'n_day4':   g['n'].mean().round(0),
    'dens':     g['dens'].mean().round(2),
    'dens_sd':  g['dens'].std().round(2),
    'cov':      g['cov'].mean().round(3),
    'd4_um':    g['d4'].mean().round(0),
    'd4_sd':    g['d4'].std().round(0),
    'd_cv':     g['d_cv'].mean().round(3),
    'log2sd':   g['log2sd'].mean().round(2),
    'dR_12':    g['g12'].mean().round(1),
    'dR_23':    g['g23'].mean().round(1),
    'dR_34':    g['g34'].mean().round(1),
})
print(out.to_string())
print()
print("  MEASURED, 8 day-4 JAPI wells:  density 0.87 to 4.11 /mm2   diameter 101 to 235 um")
print("  MEASURED, Ben's time-lapse:    density 5.6 /mm2 at 97 h    diameter 141 um")
print("  MEASURED growth increments:    +16.9, +12.5, +4.2 um/day (time-lapse);")
print("                                 BS70 +46/+33/+14, BS74 +40/+30/+21, BS77 +16/+12/0")
print("  MEASURED log2 area SD:         1.20 to 1.79 across the 8 wells")
print()
from scipy import stats
rho, p = stats.spearmanr(df.rate, df.dens)
print(f"  rate -> density      Spearman rho = {rho:+.3f}  p = {p:.2g}")
rho2, p2 = stats.spearmanr(df.dens, df.d4)
print(f"  density -> diameter  Spearman rho = {rho2:+.3f}  p = {p2:.2g}   "
      f"(across the 8 real wells this was +0.02, p = 0.96)")
if 'amp1' in df:
    print()
    print("  spatial order along the sweep:")
    print(df.groupby('rate')[['amp1', 'g030', 'S_band']].mean().round(3).to_string())
    print("  in vitro day-4 JAPI, pooled: amp1 0.140, g(0.3a) 0.301, S_band 0.60")
np.save(os.path.join(ROOT, 'results', 'nucleation_sweep', 'gr_by_rate.npy'),
        np.array([gr_by_rate], dtype=object), allow_pickle=True)
print()
print("wrote sweep_per_run.csv and gr_by_rate.npy")
