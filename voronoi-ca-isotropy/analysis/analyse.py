#!/usr/bin/env python3
"""Analyse front-propagation timings and produce the result figures.

Reads results/{lattice,voronoi}-WxH.csv written by the Lisp experiment, prints
the anisotropy summary, and writes two figures into results/.

Usage:  python3 analysis/analyse.py [results_dir]
"""

from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Categorical slots 1 and 2 of a CVD-validated palette.
# Worst all-pairs CVD dE 24.7, normal-vision dE 33.6, both >= 3:1 on a light surface.
COLOURS = {"lattice": "#2a78d6", "voronoi": "#eb6834"}
LABELS = {"lattice": "Square lattice", "voronoi": "Voronoi tessellation"}

INK = "#0b0b0b"
INK_2 = "#52514e"
GRID = "#d8d7d2"
SURFACE = "#fcfcfb"

AXIAL = ["N", "E", "S", "W"]
DIAGONAL = ["NW", "NE", "SW", "SE"]

# compass bearing -> matplotlib polar angle (radians, 0 = East, anticlockwise)
ANGLE = {
    "E": 0, "NE": 45, "N": 90, "NW": 135,
    "W": 180, "SW": 225, "S": 270, "SE": 315,
}


def ring_radii(width: int) -> tuple[int, int]:
    """Probe-ring radii in lattice cells, matching measurement.lisp."""
    end_x = width - 1
    return end_x // 4, (end_x // 10) * 4


def load(path: Path) -> tuple[list[dict], int]:
    rows = list(csv.DictReader(path.open()))
    m = re.search(r"-(\d+)x(\d+)\.csv$", path.name)
    width = int(m.group(1)) if m else 200
    return rows, width


def column(rows, name) -> np.ndarray:
    return np.array([float(r[name]) for r in rows])


def check_reached(rows) -> None:
    """Abort loudly if any probe was never reached; a 0/-1 would silently
    corrupt every mean downstream."""
    bad = [k for r in rows for k, v in r.items()
           if k not in ("replicate", "seed", "regions") and float(v) < 0]
    if bad:
        sys.exit("ERROR: front never reached probe(s) "
                 f"{sorted(set(bad))} — increase the step budget.")


def summarise(rows, ring_suffix) -> dict:
    """Mean axial / diagonal arrival time and the anisotropy index for one ring."""
    ax = np.concatenate([column(rows, d + ring_suffix) for d in AXIAL])
    dg = np.concatenate([column(rows, d + ring_suffix) for d in DIAGONAL])
    allv = np.concatenate([ax, dg])
    return {
        "axial": ax.mean(),
        "diagonal": dg.mean(),
        "anisotropy": abs(ax.mean() - dg.mean()) / allv.mean() * 100,
        "cv": allv.std() / allv.mean() * 100,
        "n": len(rows),
    }


def direction_means(rows, ring_suffix) -> dict[str, float]:
    return {d: column(rows, d + ring_suffix).mean() for d in AXIAL + DIAGONAL}


def style_axes(ax):
    ax.set_facecolor(SURFACE)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(GRID)
    ax.tick_params(colors=INK_2, length=3)
    ax.grid(color=GRID, linewidth=0.7, alpha=0.9)
    ax.set_axisbelow(True)


def figure_main(data, width, out_path):
    """Panel A: front shape as a polar speed plot. Panel B: anisotropy index."""
    r1, r2 = ring_radii(width)
    fig = plt.figure(figsize=(11, 5.6), facecolor=SURFACE)
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.15], wspace=0.30,
                          left=0.06, right=0.97, top=0.74, bottom=0.19)

    # ---- Panel A: polar front speed -------------------------------------
    axp = fig.add_subplot(gs[0, 0], projection="polar")
    axp.set_facecolor(SURFACE)
    label_anchor = {}
    for cond, rows in data.items():
        means = direction_means(rows, "")
        # speed = distance / time; normalised so the axial directions sit at 1.0
        speeds = {d: r1 / means[d] for d in means}
        base = np.mean([speeds[d] for d in AXIAL])
        theta, radius = [], []
        for d in ["E", "NE", "N", "NW", "W", "SW", "S", "SE"]:
            theta.append(np.deg2rad(ANGLE[d]))
            radius.append(speeds[d] / base)
        theta.append(theta[0])
        radius.append(radius[0])
        axp.plot(theta, radius, color=COLOURS[cond], linewidth=2,
                 marker="o", markersize=5, label=LABELS[cond],
                 markeredgecolor=SURFACE, markeredgewidth=1.2, zorder=3)
        # anchor each direct label on a different compass point so the
        # leader lines never cross the other series
        anchor_dir = "NE" if cond == "lattice" else "N"
        label_anchor[cond] = (anchor_dir, speeds[anchor_dir] / base)
    axp.plot(np.linspace(0, 2 * np.pi, 200), np.ones(200),
             color=INK_2, linewidth=1, linestyle=(0, (4, 3)), zorder=1)
    axp.set_xticks(np.deg2rad([0, 45, 90, 135, 180, 225, 270, 315]))
    axp.set_xticklabels(["E", "NE", "N", "NW", "W", "SW", "S", "SE"],
                        color=INK_2, fontsize=9)
    axp.set_ylim(0, 1.75)
    axp.set_yticks([0.5, 1.0, 1.5])
    axp.set_yticklabels(["0.5", "1.0", "1.5"], color=INK_2, fontsize=8)
    axp.set_rlabel_position(249)
    axp.grid(color=GRID, linewidth=0.7)
    axp.spines["polar"].set_color(GRID)
    for lbl in axp.get_yticklabels():
        lbl.set_bbox(dict(facecolor=SURFACE, edgecolor="none", pad=1.2))
    # direct labels, so identity never rests on colour alone
    offsets = {"lattice": (58, 0.40), "voronoi": (104, 0.34)}
    for cond, (anchor_dir, r) in label_anchor.items():
        deg, dr = offsets[cond]
        axp.annotate(LABELS[cond], xy=(np.deg2rad(ANGLE[anchor_dir]), r),
                     xytext=(np.deg2rad(deg), r + dr),
                     color=COLOURS[cond], fontsize=8.5, ha="center",
                     arrowprops=dict(arrowstyle="-", color=COLOURS[cond],
                                     linewidth=0.8, shrinkA=0, shrinkB=3))
    axp.annotate("perfectly\nisotropic", xy=(np.deg2rad(197), 1.0),
                 xytext=(np.deg2rad(192), 0.42), color=INK_2, fontsize=7.5,
                 ha="center", va="center",
                 arrowprops=dict(arrowstyle="-", color=INK_2,
                                 linewidth=0.7, shrinkA=0, shrinkB=2))
    axp.text(0.5, -0.15, "Front speed by direction (ring 1, normalised to axial)",
             transform=axp.transAxes, ha="center", color=INK, fontsize=10)

    # ---- Panel B: anisotropy index --------------------------------------
    axb = fig.add_subplot(gs[0, 1])
    style_axes(axb)
    axb.grid(axis="x", visible=False)
    rings = [("", f"Ring 1\n(r = {r1} cells)"), ("2", f"Ring 2\n(r = {r2} cells)")]
    x = np.arange(len(rings))
    bar_w = 0.34
    for k, (cond, rows) in enumerate(data.items()):
        vals = [summarise(rows, suf)["anisotropy"] for suf, _ in rings]
        pos = x + (k - 0.5) * (bar_w + 0.02)
        bars = axb.bar(pos, vals, bar_w, color=COLOURS[cond],
                       label=LABELS[cond], zorder=3)
        for b, v in zip(bars, vals):
            axb.annotate(f"{v:.1f}%", (b.get_x() + b.get_width() / 2, v),
                         xytext=(0, 4), textcoords="offset points",
                         ha="center", color=INK, fontsize=9.5)
    axb.set_xticks(x)
    axb.set_xticklabels([lbl for _, lbl in rings], color=INK_2, fontsize=9.5)
    axb.set_ylabel("Anisotropy index  (%)", color=INK_2, fontsize=9.5)
    axb.set_ylim(0, max(1.0, axb.get_ylim()[1] * 1.22))
    axb.text(0.0, -0.205, "Directional bias: |axial − diagonal| as % of mean arrival time",
             transform=axb.transAxes, ha="left", color=INK, fontsize=10)

    handles, labels = axb.get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, fontsize=9.5, labelcolor=INK_2,
               loc="upper left", bbox_to_anchor=(0.055, 0.855), ncols=2,
               handlelength=1.1, handleheight=1.1, columnspacing=1.6)

    fig.suptitle("A Voronoi tessellation removes the lattice's directional bias",
                 color=INK, fontsize=13.5, x=0.06, ha="left", y=0.965)
    fig.text(0.06, 0.905,
             f"{width}×{width} lattice · {summarise(list(data.values())[0], '')['n']} "
             f"replicates per condition · all probes equidistant from the source",
             color=INK_2, fontsize=9, ha="left")
    fig.savefig(out_path, dpi=170, facecolor=SURFACE)
    plt.close(fig)


def figure_tessellations(results_dir, width, out_path):
    """Side-by-side view of the two tessellations, with the front highlighted."""
    imgs = []
    for cond in ("lattice", "voronoi"):
        p = results_dir / f"{cond}-{width}x{width}.ppm"
        if p.exists():
            imgs.append((cond, plt.imread(p)))
    if not imgs:
        return False
    fig, axes = plt.subplots(1, len(imgs), figsize=(4.6 * len(imgs), 5.0),
                             facecolor=SURFACE)
    if len(imgs) == 1:
        axes = [axes]
    for ax, (cond, img) in zip(axes, imgs):
        ax.imshow(img)
        ax.set_title(LABELS[cond], color=INK, fontsize=11, pad=8)
        ax.set_xticks([])
        ax.set_yticks([])
        for s in ax.spines.values():
            s.set_color(GRID)
    fig.suptitle("The swept region, at the moment the front reaches ring 1",
                 color=INK, fontsize=12.5, y=0.975)
    fig.text(0.5, 0.915,
             "dark = already swept   ·   red = active front   ·   "
             "pale mosaic = untouched regions",
             color=INK_2, fontsize=9, ha="center")
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    fig.savefig(out_path, dpi=150, facecolor=SURFACE)
    plt.close(fig)
    return True


def main():
    results_dir = Path(sys.argv[1] if len(sys.argv) > 1 else "results")
    files = sorted(results_dir.glob("*-*x*.csv"))
    if not files:
        sys.exit(f"No result CSVs in {results_dir}/ — run the experiment first.")

    data, width = {}, None
    for f in files:
        cond = f.name.split("-")[0]
        rows, width = load(f)
        data[cond] = rows
    data = {k: data[k] for k in ("lattice", "voronoi") if k in data}
    for rows in data.values():
        check_reached(rows)

    print(f"\n{width}×{width} lattice\n")
    hdr = f"{'condition':<10}{'ring':<7}{'regions':>9}{'axial':>9}{'diagonal':>10}{'anisotropy':>12}{'CV':>8}"
    print(hdr)
    print("-" * len(hdr))
    for cond, rows in data.items():
        regions = int(column(rows, "regions").mean())
        for suf, name in (("", "1"), ("2", "2")):
            s = summarise(rows, suf)
            print(f"{cond:<10}{name:<7}{regions:>9}{s['axial']:>9.2f}"
                  f"{s['diagonal']:>10.2f}{s['anisotropy']:>11.2f}%{s['cv']:>7.1f}%")
    print()
    if "lattice" in data and "voronoi" in data:
        for suf, name in (("", "1"), ("2", "2")):
            l = summarise(data["lattice"], suf)["anisotropy"]
            v = summarise(data["voronoi"], suf)["anisotropy"]
            factor = f"{l / v:.1f}×" if v > 0.01 else "n/a"
            print(f"  ring {name}: anisotropy {l:.2f}% → {v:.2f}%  "
                  f"(reduced {l - v:.2f} pp, {factor})")

    results_dir.mkdir(exist_ok=True)
    main_fig = results_dir / "anisotropy.png"
    figure_main(data, width, main_fig)
    print(f"\nwrote {main_fig}")
    tess_fig = results_dir / "tessellations.png"
    if figure_tessellations(results_dir, width, tess_fig):
        print(f"wrote {tess_fig}")


if __name__ == "__main__":
    main()
