"""Summary figure of the ladder batch (REPORT section 13): analysis/summaries/ladder_summary.png. Usage: python analysis/ladder_summary_fig.py"""

import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FixedFormatter, NullLocator
order = ["bi13.5","bi14.1","bi15","bi20","bi25","bi30","bi50"]; bi = [13.5,14.1,15,20,25,30,50]; xs = np.arange(len(bi))
conds = ["seed9","seed36","seed144","random","nucl"]; labels = {"seed9":"9 seeds","seed36":"36 seeds","seed144":"144 seeds","random":"random U(0,2)","nucl":"nucleation"}
col = dict(zip(conds, ["#2a78d6","#eb6834","#1baf7a","#eda100","#e87ba4"]))
f = pd.read_csv("analysis/summaries/ladder_frames.csv"); f = f[(f.channel=="a") & np.isclose(f.t, 600.01)]
x = pd.read_csv("analysis/summaries/ladder_fft.csv"); x = x[(x.channel=="a") & (x.t==600)]
def agg(d, cond, colname, fn="mean"):
    s = d[d.condition==cond].groupby("set")[colname]
    return np.array([getattr(s, fn)().get(k, np.nan) for k in order])
fig, ax = plt.subplots(2, 2, figsize=(10.5, 8.2)); fig.patch.set_facecolor("#fcfcfb")
for a in ax.ravel():
    a.set_facecolor("#fcfcfb"); a.grid(True, color="#e6e5e1", linewidth=0.6); a.set_axisbelow(True)
    for sp in ("top", "right"): a.spines[sp].set_visible(False)
    a.spines["left"].set_color("#c3c2b7"); a.spines["bottom"].set_color("#c3c2b7"); a.tick_params(colors="#52514e", labelsize=8.5)
    a.set_xticks(xs); a.set_xticklabels([f"{b:g}" for b in bi]); a.set_xlabel("b_i (ordinal axis)", fontsize=9, color="#52514e")
def lines(a, ycol, src, yl, logy=False, title="", conds_=conds, lone=None, legend_kw=None):
    for c in conds_:
        y = agg(src, c, ycol); lo = agg(src, c, ycol, "min"); hi = agg(src, c, ycol, "max")
        a.plot(xs, y, "-", color=col[c], linewidth=2, marker="o", markersize=5, markerfacecolor="#fcfcfb", markeredgewidth=1.6, label=labels[c])
        a.fill_between(xs, lo, hi, color=col[c], alpha=0.12, linewidth=0)
    if lone is not None:
        a.plot(xs, lone, "--", color="#0b0b0b", linewidth=1.2, label="lone seed")
    if logy:
        a.set_yscale("log"); a.yaxis.set_minor_locator(NullLocator())
    a.set_ylabel(yl, fontsize=9, color="#52514e"); a.set_title(title, fontsize=10, loc="left", color="#0b0b0b")
    a.legend(fontsize=8, frameon=False, labelcolor="#52514e", **(legend_kw or {}))
lone = agg(f, "seed1", "area_real_median")
lines(ax[0,0], "area_real_median", f, "median real-domain area (cells)", logy=True, title="a  Domain size at t = 600 (bands: replicate range)", lone=lone, legend_kw=dict(loc="lower left"))
ax[0,0].yaxis.set_major_locator(FixedLocator([10, 20, 50, 100, 200, 500, 1000])); ax[0,0].yaxis.set_major_formatter(FixedFormatter(["10","20","50","100","200","500","1000"]))
ratio = f.assign(ratio=f.nn_cv_border / f.csr_cv_border)
lines(ax[0,1], "ratio", ratio, "border-corrected NN CV / same-n CSR reference", title="b  Regularity relative to CSR (1 = random; 9-seed fields omitted)", conds_=["seed36","seed144","random","nucl"], legend_kw=dict(loc="upper left"))
ax[0,1].axhline(1.0, color="#c3c2b7", linewidth=1); ax[0,1].set_ylim(0, 1.1)
lines(ax[1,0], "fft_peak_wavelength", x, "spectral peak wavelength (cells)", title="c  Hex-correct power-spectrum peak at t = 600", conds_=["seed144","random","nucl"], legend_kw=dict(loc="upper center"))
ax[1,0].axhline(7.75, color="#c3c2b7", linewidth=1); ax[1,0].annotate("inhibitor range sqrt(1.5 D_i/gamma) = 7.75", (0, 7.75), xytext=(2, 3), textcoords="offset points", fontsize=8, color="#52514e"); ax[1,0].set_ylim(0, 33)
dens = f.assign(per1e4=f.blobs_real / 4.0)
lines(ax[1,1], "per1e4", dens, "real domains per 10,000 cells", logy=True, title="d  Domain density at t = 600", legend_kw=dict(loc="center left", bbox_to_anchor=(0.0, 0.42)))
ax[1,1].yaxis.set_major_locator(FixedLocator([2, 5, 10, 20, 50])); ax[1,1].yaxis.set_major_formatter(FixedFormatter(["2","5","10","20","50"]))
fig.suptitle("Ladder batch at n 10/4, b_a 5, D_i 20, gamma 0.5: 200x200, six initial conditions, t = 600 (3 replicates per field condition)", fontsize=10, color="#0b0b0b")
fig.tight_layout(rect=(0, 0, 1, 0.96)); fig.savefig("analysis/summaries/ladder_summary.png", dpi=160, facecolor=fig.get_facecolor()); print("written")
