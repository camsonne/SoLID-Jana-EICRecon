#!/usr/bin/env python3
"""Plot the DDVCS missing-mass luminosity scan produced by run_ddvcs_lumi_scan.sh.

    plot_ddvcs_results.py RESULTS_DIR

Reads RESULTS_DIR/L<lumi>_<readout>_summary.json and *_events.csv, writes PNG figures and
RESULTS.md (table view of every number in the plots) into RESULTS_DIR.
"""
import csv
import glob
import json
import math
import os
import re
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

SURFACE = "#fcfcfb"
TEXT = "#0b0b0b"
TEXT2 = "#52514e"
GRID = "#e6e5e1"
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]  # fixed categorical order by luminosity
MP2 = 0.938272088 ** 2


def style(ax):
    ax.set_facecolor(SURFACE)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.tick_params(colors=TEXT2, labelsize=9)
    ax.yaxis.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)


def fmt_lumi(L):
    e = int(math.floor(math.log10(L)))
    m = L / 10 ** e
    return (f"{m:.1f}" if abs(m - round(m)) > 0.05 else f"{m:.0f}") + rf"$\times10^{{{e}}}$"


def load(resdir):
    runs = []
    for f in sorted(glob.glob(os.path.join(resdir, "*_summary.json"))):
        j = json.load(open(f))
        m = re.match(r"L([0-9.e+]+)_(\w+)_summary\.json", os.path.basename(f))
        if not m:
            continue
        L = float(m.group(1))
        ro = m.group(2)
        csvf = f.replace("_summary.json", "_events.csv")
        mm2 = []
        if os.path.exists(csvf):
            for r in csv.DictReader(open(csvf)):
                if r["found"] == "1":
                    mm2.append(float(r["mm2"]))
        occ = {}
        logf = f.replace("_summary.json", ".log")
        if os.path.exists(logf):
            for line in open(logf):
                mo = re.search(r"GEM(\d) @Rin: rate=([0-9.e+-]+) kHz/mm\^2, hits=([0-9.e+-]+)/cm\^2, occupancy=([0-9.e+-]+), survival=([0-9.e+-]+)", line)
                if mo:
                    occ[int(mo.group(1))] = tuple(float(mo.group(i)) for i in range(2, 6))
        runs.append(dict(L=L, readout=ro, summary=j, mm2=mm2, occ=occ))
    runs.sort(key=lambda r: (r["readout"], r["L"]))
    return runs


def plot_mm2(runs, resdir):
    readouts = sorted({r["readout"] for r in runs}, reverse=True)  # strip first
    lumis = sorted({r["L"] for r in runs})
    fig, axes = plt.subplots(1, len(readouts), figsize=(5.2 * len(readouts), 4.2), squeeze=False, facecolor=SURFACE)
    bins = [0.4 + 0.01 * i for i in range(101)]
    for ax, ro in zip(axes[0], readouts):
        style(ax)
        for r in [x for x in runs if x["readout"] == ro]:
            ci = lumis.index(r["L"]) % len(SERIES)
            s = r["summary"]["mm2"]
            if not r["mm2"]:
                continue
            label = f"L = {fmt_lumi(r['L'])} cm$^{{-2}}$s$^{{-1}}$: $\\sigma$ = {s['gauss_sigma']*1e3:.0f} MeV$^2$"
            ax.hist(r["mm2"], bins=bins, histtype="step", linewidth=1.8, color=SERIES[ci], density=True, label=label)
        ax.axvline(MP2, color=TEXT2, linewidth=1, linestyle=":")
        ax.text(MP2 + 0.01, ax.get_ylim()[1] * 0.95, "$M_p^2$", color=TEXT2, fontsize=9, va="top")
        ax.set_xlabel("Missing mass squared  $MM^2 = (k + P - k' - \\mu^+ - \\mu^-)^2$  [GeV$^2$]", color=TEXT, fontsize=9.5)
        ax.set_ylabel("Probability density", color=TEXT, fontsize=9.5)
        ax.set_title(f"GEM {ro} readout", color=TEXT, fontsize=11, loc="left")
        ax.legend(frameon=False, fontsize=8.5, labelcolor=TEXT)
    fig.suptitle("SoLID DDVCS  e p → e' p (μ⁺μ⁻): missing-mass distribution of the undetected proton", color=TEXT, fontsize=11, x=0.02, ha="left")
    fig.tight_layout()
    fig.savefig(os.path.join(resdir, "mm2_distributions.png"), dpi=150, facecolor=SURFACE)
    plt.close(fig)


def plot_vs_lumi(runs, resdir):
    readouts = sorted({r["readout"] for r in runs}, reverse=True)
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), facecolor=SURFACE)
    styles = {"strip": dict(linestyle="-", marker="o"), "pixel": dict(linestyle="--", marker="s")}
    for i, ro in enumerate(readouts):
        rs = [x for x in runs if x["readout"] == ro]
        L = [x["L"] for x in rs]
        col = SERIES[i % len(SERIES)]
        st = styles.get(ro, dict(linestyle="-", marker="o"))
        sg = [x["summary"]["mm2"]["gauss_sigma"] * 1e3 for x in rs]
        s68 = [x["summary"]["mm2"]["sigma68"] * 1e3 for x in rs]
        axes[0].plot(L, sg, color=col, linewidth=2, markersize=7, label=f"{ro}: Gaussian core", **st)
        axes[0].plot(L, s68, color=col, linewidth=1, markersize=7, markerfacecolor=SURFACE, alpha=0.8, label=f"{ro}: $\\sigma_{{68}}$", **st)
        for x, y in zip(L, sg):
            axes[0].annotate(f"{y:.0f}", (x, y), textcoords="offset points", xytext=(0, 7), ha="center", fontsize=8, color=TEXT2)
        eff = [x["summary"]["reco_efficiency_in_acceptance"] for x in rs]
        axes[1].plot(L, eff, color=col, linewidth=2, markersize=7, label=ro, **st)
        for x, y in zip(L, eff):
            axes[1].annotate(f"{y:.2f}", (x, y), textcoords="offset points", xytext=(0, 7), ha="center", fontsize=8, color=TEXT2)
        fb = [x["summary"]["frac_tracks_with_bkg_hit"] for x in rs]
        axes[2].plot(L, fb, color=col, linewidth=2, markersize=7, label=ro, **st)
    for ax in axes:
        style(ax)
        ax.set_xscale("log")
        ax.set_xlabel("Luminosity [cm$^{-2}$ s$^{-1}$]", color=TEXT, fontsize=9.5)
        ax.legend(frameon=False, fontsize=8.5, labelcolor=TEXT)
    axes[0].set_ylabel("$\\sigma(MM^2)$  [MeV$^2$]", color=TEXT, fontsize=9.5)
    axes[0].set_title("Missing-mass resolution", color=TEXT, fontsize=11, loc="left")
    axes[0].set_ylim(bottom=0)
    axes[1].set_ylabel("Fraction of accepted events with e', μ⁺, μ⁻ reconstructed", color=TEXT, fontsize=9)
    axes[1].set_title("Reconstruction efficiency", color=TEXT, fontsize=11, loc="left")
    axes[1].set_ylim(0, 1.05)
    axes[2].set_ylabel("Fraction of tracks containing a background hit", color=TEXT, fontsize=9)
    axes[2].set_title("Background contamination of tracks", color=TEXT, fontsize=11, loc="left")
    axes[2].set_ylim(bottom=0)
    fig.tight_layout()
    fig.savefig(os.path.join(resdir, "resolution_vs_luminosity.png"), dpi=150, facecolor=SURFACE)
    plt.close(fig)


def write_table(runs, resdir):
    lines = ["# DDVCS missing-mass resolution vs luminosity", "",
             "Generated by `scripts/plot_ddvcs_results.py` from the fast-simulation + JANA2 reconstruction chain.", "",
             "| L [cm⁻² s⁻¹] | GEM readout | events | acceptance | ε(reco \\| acc) | σ(MM²) core [MeV²] | σ₆₈(MM²) [MeV²] | ⟨MM²⟩ [GeV²] | σ(M_X) core [MeV] | σ(p)/p μ | tracks w/ bkg hit | fake tracks | GEM2 occupancy @Rin |",
             "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in runs:
        s = r["summary"]
        occ = r["occ"].get(2, (0, 0, 0, 0))[2]
        lines.append(
            f"| {r['L']:.2g} | {r['readout']} | {s['n_events']} | {s['acceptance']:.3f} | {s['reco_efficiency_in_acceptance']:.3f} | "
            f"{s['mm2']['gauss_sigma']*1e3:.0f} | {s['mm2']['sigma68']*1e3:.0f} | {s['mm2']['gauss_mean']:.3f} | {s['mm']['gauss_sigma']*1e3:.1f} | "
            f"{s['dp_over_p_muon']['gauss_sigma']*100:.2f}% | {s['frac_tracks_with_bkg_hit']:.3f} | {s['frac_fake_tracks']:.4f} | {occ:.3f} |")
    lines += ["", "## GEM background at the inner radius of each plane", "",
              "| L | readout | plane | rate [kHz/mm²] | hits in 100 ns [1/cm²] | channel occupancy | hit survival |", "|---|---|---|---|---|---|---|"]
    for r in runs:
        for pl in sorted(r["occ"]):
            rate, hits, occ, surv = r["occ"][pl]
            lines.append(f"| {r['L']:.2g} | {r['readout']} | {pl} | {rate:.1f} | {hits:.3f} | {occ:.3f} | {surv:.3f} |")
    lines += ["", "![MM2 distributions](mm2_distributions.png)", "", "![Resolution vs luminosity](resolution_vs_luminosity.png)", ""]
    open(os.path.join(resdir, "RESULTS.md"), "w").write("\n".join(lines))


def main():
    resdir = sys.argv[1] if len(sys.argv) > 1 else "results"
    runs = load(resdir)
    if not runs:
        print("no *_summary.json found in", resdir)
        return 1
    plot_mm2(runs, resdir)
    plot_vs_lumi(runs, resdir)
    write_table(runs, resdir)
    print("wrote", os.path.join(resdir, "RESULTS.md"), "and figures")
    return 0


if __name__ == "__main__":
    sys.exit(main())
