#!/usr/bin/env python3
"""Compare two luminosity scans (e.g. placeholder vs analytic background table).

    compare_scans.py results results_analytic [--labels placeholder analytic] [-o results_analytic/COMPARISON.md]
"""
import argparse
import glob
import json
import os
import re


def load(resdir):
    runs = {}
    for f in sorted(glob.glob(os.path.join(resdir, "*_summary.json"))):
        m = re.match(r"L([0-9.e+]+)_(\w+)_summary\.json", os.path.basename(f))
        if not m:
            continue
        runs[(float(m.group(1)), m.group(2))] = json.load(open(f))
    return runs


def fmt(x, nd=3):
    return "—" if x is None else f"{x:.{nd}f}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dir_a")
    ap.add_argument("dir_b")
    ap.add_argument("--labels", nargs=2, default=["placeholder", "analytic"])
    ap.add_argument("-o", "--output", default=None)
    args = ap.parse_args()
    A, B = load(args.dir_a), load(args.dir_b)
    la, lb = args.labels
    keys = sorted(set(A) | set(B), key=lambda k: (k[1] != "strip", k[0]))
    lines = [f"# Luminosity scan: {la} vs {lb} background table", "",
             f"`{la}` = `{args.dir_a}`, `{lb}` = `{args.dir_b}`. Efficiency is the fraction of events with e′, μ⁺, μ⁻ "
             "all inside the geometric acceptance that are fully reconstructed; σ(MM²) is the Gaussian-core width.", "",
             f"| L [cm⁻²s⁻¹] | readout | events | ε {la} | ε {lb} | σ(MM²) {la} [MeV²] | σ(MM²) {lb} [MeV²] | "
             f"tracks w/ bkg hit {la} | {lb} | fakes {la} | {lb} |",
             "|---|---|---|---|---|---|---|---|---|---|---|"]
    for k in keys:
        a, b = A.get(k), B.get(k)
        def g(s, *path, scale=1.0, nd=3):
            if s is None:
                return None
            v = s
            for p in path:
                v = v.get(p) if isinstance(v, dict) else None
                if v is None:
                    return None
            return v * scale
        def eff(s):
            if s is None:
                return None
            return s.get("efficiency_vs_geometric_acceptance", s.get("reco_efficiency_in_acceptance"))
        def sig(s):
            v = g(s, "mm2", "gauss_sigma", scale=1e3, nd=0)
            n = g(s, "mm2", "n")
            return None if not n else v
        nev = (a or b)["n_events"]
        lines.append(
            f"| {k[0]:.2g} | {k[1]} | {nev} | {fmt(eff(a))} | {fmt(eff(b))} | {fmt(sig(a), 0)} | {fmt(sig(b), 0)} | "
            f"{fmt(g(a, 'frac_tracks_with_bkg_hit'))} | {fmt(g(b, 'frac_tracks_with_bkg_hit'))} | "
            f"{fmt(g(a, 'frac_fake_tracks'), 4)} | {fmt(g(b, 'frac_fake_tracks'), 4)} |")
    text = "\n".join(lines) + "\n"
    if args.output:
        open(args.output, "w").write(text)
        print("wrote", args.output)
    print(text)


if __name__ == "__main__":
    main()
