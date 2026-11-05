"""Redraw the three figures from data/*.csv. Usage: python figures.py  (needs matplotlib, numpy, pandas)."""
import os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER, INK, ORANGE = "#f4efe4", "#222222", "#d9622b"
GREYS = {"dark": "#555555", "mid": "#8a8a8a", "light": "#b5b0a4"}
XLAB, YLAB = "steps", "pull toward light (mean PT)"
plt.rcParams.update({
    "figure.facecolor": PAPER, "axes.facecolor": PAPER, "savefig.facecolor": PAPER,
    "text.color": INK, "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK,
    "axes.edgecolor": INK, "axes.linewidth": 0.6, "font.family": "serif", "font.size": 15,
    "axes.spines.top": False, "axes.spines.right": False,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
})

def load(name):
    return pd.read_csv(os.path.join(HERE, "data", name))

def new_fig(ncols=1):
    fig, axes = plt.subplots(1, ncols, figsize=(16, 10), dpi=100, squeeze=False)
    return fig, axes[0]

def draw(ax, d, run, color, lw, col="pt_mean", x0=0):
    r = d[(d.run == run) & (d.step > x0)].dropna(subset=[col])
    ax.plot(r.step, r[col], color=color, lw=lw, solid_capstyle="round")
    return r

def label(ax, text, x, y, color, dy=0):
    ax.annotate(text, (x, y), xytext=(8, dy), textcoords="offset points", color=color,
                va="center", ha="left", fontsize=15, annotation_clip=False)

def groups(d, arms):
    """arms: list of (prefix, color, lw, label dy in points). Seeds are drawn thin, the arm named once at the line end."""
    for prefix, color, lw, dy in arms:
        last = []
        for run in sorted(r for r in d.run.unique() if r.startswith(prefix) and r[len(prefix):].isdigit()):
            r = draw(ax_cur, d, run, color, lw)
            last.append((r.step.iloc[-1], r.pt_mean.iloc[-1]))
        x = max(p[0] for p in last); y = float(np.mean([p[1] for p in last]))
        label(ax_cur, prefix, x, y, color, dy)

def finish(fig, ax, title, path):
    ax.set_xlabel(XLAB); ax.set_ylabel(YLAB)
    ax.set_ylim(0, 1.0)
    ax.grid(axis="y", color=INK, alpha=0.12, lw=0.5)
    fig.text(0.07, 0.955, title, fontsize=19, color=INK, ha="left")
    fig.subplots_adjust(left=0.08, right=0.9, top=0.9, bottom=0.1, wspace=0.3)
    fig.savefig(os.path.join(HERE, path), dpi=100)
    plt.close(fig)

# fig 1: the moving sky, three arms of three seeds, sampled every 150,000 steps
d = load("moving_sky.csv")
fig, (ax_cur,) = new_fig()
groups(d, [("N", GREYS["dark"], 1.2, 0), ("E", GREYS["light"], 1.2, 0), ("H", ORANGE, 1.6, 0)])
ax_cur.set_xlim(0, 3.0e6); ax_cur.xaxis.set_major_formatter(lambda v, p: f"{v/1e6:g}M")
finish(fig, ax_cur, "The moving sky: with hunters (H), without (N), light changing over the whole glass (E)", "fig1_moving_sky.png")

# fig 2: the four fixed-sky arms, three seeds each, log steps so the early fall shows
d = load("fixed_sky.csv")
fig, (ax_cur,) = new_fig()
groups(d, [("UN", GREYS["light"], 1.0, 14), ("UH", GREYS["mid"], 1.0, -14), ("GN", GREYS["dark"], 1.3, 0), ("GH", ORANGE, 1.6, 0)])
ax_cur.set_xscale("log"); ax_cur.set_xlim(3e3, 1.5e6)
ax_cur.xaxis.set_major_formatter(lambda v, p: f"{v:,.0f}")
finish(fig, ax_cur, "A fixed sky: G is the gradient, U is even light; H has hunters, N has none", "fig2_fixed_sky.png")

# fig 3: the density arm beside the hunter and hunterless arms, and where the bodies sit
d = load("density.csv")
fig, axes = new_fig(2)
ax_cur = axes[0]
groups(d, [("GN", GREYS["dark"], 1.3, 8), ("GD", GREYS["mid"], 1.3, -8), ("GH", ORANGE, 1.6, 0)])
ax_cur.set_xscale("log"); ax_cur.set_xlim(3e3, 1.5e6)
ax_cur.xaxis.set_major_formatter(lambda v, p: f"{v:,.0f}")
ax_cur.set_xlabel(XLAB); ax_cur.set_ylabel(YLAB); ax_cur.set_ylim(0, 1.0)
ax_cur.grid(axis="y", color=INK, alpha=0.12, lw=0.5)
ax_cur.set_title("pull toward light", loc="left", fontsize=15)
ax_cur = axes[1]
for prefix, color, lw, dy in [("GN", GREYS["dark"], 1.3, 8), ("GD", GREYS["mid"], 1.3, -8), ("GH", ORANGE, 1.6, 0)]:
    ends = []
    for run in sorted(r for r in d.run.unique() if r.startswith(prefix)):
        r = draw(ax_cur, d, run, color, lw, col="top5_share")
        ends.append((r.step.iloc[-1], r.top5_share.iloc[-1]))
    label(ax_cur, prefix, max(e[0] for e in ends), float(np.mean([e[1] for e in ends])), color, dy)
ax_cur.set_xscale("log"); ax_cur.set_xlim(1.4e5, 1.5e6)
ax_cur.set_xticks([2e5, 5e5, 1e6]); ax_cur.xaxis.set_major_formatter(lambda v, p: f"{v:,.0f}")
ax_cur.xaxis.set_minor_formatter(lambda v, p: "")
ax_cur.set_xlabel(XLAB); ax_cur.set_ylabel("share of bodies in the five brightest columns"); ax_cur.set_ylim(0, 1.0)
ax_cur.grid(axis="y", color=INK, alpha=0.12, lw=0.5)
ax_cur.set_title("where the bodies sit", loc="left", fontsize=15)
fig.text(0.07, 0.955, "Density arm GD: three times the upkeep, no hunters, the most crowded at the light", fontsize=19, ha="left")
fig.subplots_adjust(left=0.07, right=0.94, top=0.88, bottom=0.1, wspace=0.28)
fig.savefig(os.path.join(HERE, "fig3_density.png"), dpi=100)
