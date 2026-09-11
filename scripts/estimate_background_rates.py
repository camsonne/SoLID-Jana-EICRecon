#!/usr/bin/env python3
"""Analytic estimate of the GEM background rate densities for the SoLID DDVCS configuration.

The estimate folds textbook electromagnetic and hadronic cross sections (in the form used by the
Geant4 Physics Reference Manual and the PDG) with the actual target, beam current, solenoid field
at the target and GEM geometry. It is a *direct-source floor*: it contains the particles that
travel from the target region straight to the GEM planes. What it cannot contain is the albedo
from the downstream beam pipe, collimators and endcap, into which the multi-THz Moller flux and the
~500 W of bremsstrahlung photons are dumped; that term needs a Geant4 beam-on-target simulation
(see docs/background_estimate.md).

Sources included
  1. Moller electrons e- e- -> e- e- on target electrons  [G4 PRM: Moller-Bhabha ionisation]
     Direct hits only if the electron is not confined by the solenoid or reflected by the
     magnetic mirror (helix integrated numerically in the on-axis field).
  2. Bremsstrahlung photons radiated by those Moller electrons while leaving the LH2 and the Al
     cell wall  [G4 PRM: Seltzer-Berger / Bethe-Heitler, complete screening]
  3. Compton scattering (Klein-Nishina) of the primary-beam bremsstrahlung photons inside the target
     [G4 PRM: Klein-Nishina]
  4. Hadronic photo/electro-production (equivalent-photon flux x sigma_gamma-p total, PDG) with an
     assumed forward-peaked angular distribution of the secondaries.
Photon hits are weighted with the interaction probability in one GEM plane (Compton + photoelectric
in the Cu + pair production) and a 50% probability that an interaction gives a signal.

Usage
  scripts/estimate_background_rates.py [--luminosity 1.2e37] [--target-length 15] [--window-ns 100]
      [--flat-albedo-kHz-mm2 0] [--config-out config/gem_background_rates_analytic.txt]
"""
import argparse
import math
import sys

import numpy as np

_trapz = getattr(np, "trapezoid", None) or np.trapz  # numpy 2.x renamed trapz

# ----------------------------------------------------------------------------------------------
# constants (MeV, cm, s)
ALPHA = 1 / 137.036
R_E = 2.8179403e-13          # cm
M_E = 0.51099895             # MeV
N_A = 6.02214076e23
E_CHARGE = 1.602176634e-19   # C
BARN = 1e-24                 # cm^2

# ----------------------------------------------------------------------------------------------
# SoLID DDVCS geometry (identical to SoLIDGeometry_service)
TARGET_Z = -315.0            # cm
TARGET_RADIUS = 1.8872       # cm (LH2)
WALL_AL = 0.0178             # cm
LH2_RHO = 0.0708             # g/cm3
LH2_X0 = 890.4               # cm
AL_X0 = 8.897                # cm
N_E_LH2 = LH2_RHO * N_A / 1.00794   # electrons (= protons) per cm3
PLANES = [  # id, z, rin, rout
    (1, -175.0, 39.0, 135.0), (2, -150.0, 21.0, 135.0), (3, -119.0, 25.0, 135.0),
    (4, -68.0, 32.0, 135.0), (5, 5.0, 42.0, 99.0), (6, 92.0, 53.0, 122.0),
]
PLACEHOLDER = {1: 1.0, 2: 1.0, 3: 0.8, 4: 0.6, 5: 0.5, 6: 0.4}  # kHz/mm^2 at R0=40 cm, alpha=2, L=1.2e37
# GEM plane: ~0.26 g/cm^2 (Al windows, G10, Nomex, Kapton, Cu, gas), Cu mass fraction ~0.17
GEM_MASS_THICKNESS = 0.26    # g/cm^2
GEM_CU_FRACTION = 0.17
GEM_SIGNAL_PROB = 0.5        # probability that a photon interaction yields a hit

# CLEO-II solenoid, on-axis field of a thin solenoid (a = 153 cm, z in [-173.8, 173.8]) scaled to 1.5 T
COIL_A, COIL_HALF_L, B_CENTER = 153.0, 173.8, 1.5


def bz_axis(z):
    """On-axis Bz [T] of the finite solenoid (same normalisation as SoLIDField_service)."""
    def f(zz):
        return zz / math.sqrt(zz * zz + COIL_A * COIL_A)
    raw = 0.5 * (f(z + COIL_HALF_L) - f(z - COIL_HALF_L))
    raw0 = 0.5 * (f(COIL_HALF_L) - f(-COIL_HALF_L))
    return B_CENTER * raw / raw0


# ----------------------------------------------------------------------------------------------
# 1. Moller scattering
def moller_dsigma_dT(T, K):
    """dsigma/dT per target electron [cm^2/MeV]; T = recoil kinetic energy, K = beam kinetic energy."""
    gamma = 1 + K / M_E
    beta2 = 1 - 1 / gamma**2
    Tp = K - T
    return (2 * math.pi * R_E**2 * M_E / beta2) * (
        1 / T**2 + 1 / Tp**2 + ((gamma - 1) / (gamma * K))**2 - (2 * gamma - 1) / (gamma**2 * T * Tp))


def moller_T_of_theta(theta, E0):
    """Recoil kinetic energy [MeV] of a free electron at rest scattered to lab angle theta."""
    p0 = math.sqrt(E0**2 - M_E**2)
    c = math.cos(theta)
    return 2 * M_E * p0**2 * c**2 / ((E0 + M_E)**2 - p0**2 * c**2)


def moller_dN_dOmega(theta, E0, electrons_per_s):
    """Moller electrons per second per steradian at lab angle theta (per unit target areal density)."""
    K = E0 - M_E
    dth = 1e-5
    T1, T2 = moller_T_of_theta(theta - dth, E0), moller_T_of_theta(theta + dth, E0)
    dT_dtheta = abs(T2 - T1) / (2 * dth)
    T = moller_T_of_theta(theta, E0)
    return electrons_per_s * moller_dsigma_dT(T, K) * dT_dtheta / (2 * math.pi * math.sin(theta)), T


def helix_reaches(theta, T, z_plane, rin, rout, charge=-1, steps=4000):
    """Integrate a charged particle from the target centre in the on-axis field Bz(z) and report the
    radius at z_plane (None if it turns around or reflects). Field taken axial and uniform in r."""
    p = math.sqrt(T * (T + 2 * M_E)) * 1e-3   # GeV
    px, py, pz = p * math.sin(theta), 0.0, p * math.cos(theta)
    x, y, z = 0.0, 0.0, TARGET_Z
    ds = (z_plane - TARGET_Z) / steps * 1.2
    k = 0.00299792458  # GeV/(T cm)
    for _ in range(steps * 3):
        if z >= z_plane:
            return math.hypot(x, y)
        B = bz_axis(z)
        pm = math.sqrt(px * px + py * py + pz * pz)
        ux, uy, uz = px / pm, py / pm, pz / pm
        if uz <= 0:
            return None  # reflected by the mirror
        # dp/ds = q k (u x B), B along z
        dpx, dpy = charge * k * (uy * B), charge * k * (-ux * B)
        # adaptive step: fraction of the gyroradius
        pt = math.hypot(px, py)
        rg = pt / (k * max(B, 1e-3))
        h = min(ds, 0.1 * rg) if rg > 0 else ds
        x += ux * h; y += uy * h; z += uz * h
        px += dpx * h; py += dpy * h
        # renormalise |p|
        pn = math.sqrt(px * px + py * py + pz * pz)
        px *= pm / pn; py *= pm / pn; pz *= pm / pn
    return None


# ----------------------------------------------------------------------------------------------
# 2. bremsstrahlung (thin target, complete screening): photons per electron with k > kmin over t X0
def brems_yield(T, t_x0, kmin):
    if T <= kmin:
        return 0.0
    return t_x0 * (4.0 / 3.0) * (math.log(T / kmin) - (1 - kmin / T)) + t_x0 * 0.5 * (1 - kmin / T)**2 * 0


def brems_spectrum(k, E0, t_x0):
    """dN/dk [1/MeV] per electron for the primary beam over t_x0 radiation lengths."""
    y = k / E0
    return t_x0 / k * (4.0 / 3.0 - 4.0 / 3.0 * y + y * y)


# ----------------------------------------------------------------------------------------------
# 3. Compton (Klein-Nishina)
def kn_dsigma_dOmega(k, theta):
    kp = k / (1 + (k / M_E) * (1 - math.cos(theta)))
    r = kp / k
    return 0.5 * R_E**2 * r * r * (r + 1 / r - math.sin(theta)**2), kp


def kn_sigma_total(k):
    e = k / M_E
    if e < 1e-3:
        return 8 * math.pi / 3 * R_E**2
    return 2 * math.pi * R_E**2 * ((1 + e) / e**2 * (2 * (1 + e) / (1 + 2 * e) - math.log(1 + 2 * e) / e)
                                   + math.log(1 + 2 * e) / (2 * e) - (1 + 3 * e) / (1 + 2 * e)**2)


def gem_interaction_prob(k):
    """Probability that a photon of energy k [MeV] interacts in one GEM plane."""
    mu_compton = N_A * 0.5 * kn_sigma_total(k)                       # cm^2/g, Z/A = 0.5
    mu_pe_cu = 0.45 * (0.1 / k)**3 if k > 0.1 else 0.45 * (0.1 / max(k, 0.03))**2.5
    mu_pair = (7.0 / 9.0) / 25.0 * (1 - math.exp(-k / 20.0))       # X0(mixture) ~ 25 g/cm^2
    mu = mu_compton + GEM_CU_FRACTION * mu_pe_cu + mu_pair
    return 1 - math.exp(-mu * GEM_MASS_THICKNESS)


# ----------------------------------------------------------------------------------------------
# 4. hadronic photoproduction
def sigma_gamma_p(k_gev):
    """Total gamma-p cross section [cm^2] (resonance region approximated, DL fit above 2 GeV)."""
    if k_gev < 0.15:
        return 0.0
    if k_gev < 0.5:
        return 350e-6 * 1e-24
    if k_gev < 2.0:
        return 200e-6 * 1e-24
    return (0.0677 * k_gev**0.0808 + 0.129 * k_gev**-0.4525) * 1e-27


def equivalent_photon_flux(k, E0):
    """Weizsacker-Williams electron -> quasi-real photon flux dn/dk [1/MeV], Q2max ~ 1 GeV^2."""
    y = k / E0
    if y <= 0 or y >= 1:
        return 0.0
    q2min = (M_E * 1e-3)**2 * y * y / (1 - y)  # GeV^2
    lnq = math.log(1.0 / q2min)
    return ALPHA / (math.pi * k) * ((1 + (1 - y)**2) * lnq - 2 * (1 - y))


def hadron_angular(theta):
    """Assumed angular distribution of hadronic secondaries, dN/dOmega ~ 1/(theta^2 + thc^2), normalised
    to 1 over the forward hemisphere."""
    thc = math.radians(20.0)
    f = lambda th: 1.0 / (th * th + thc * thc)
    ths = np.linspace(1e-3, math.pi / 2, 2000)
    norm = _trapz([f(t) * 2 * math.pi * math.sin(t) for t in ths], ths)
    return f(theta) / norm


# ----------------------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--luminosity", type=float, default=1.2e37, help="cm^-2 s^-1")
    ap.add_argument("--beam-energy", type=float, default=11000.0, help="MeV")
    ap.add_argument("--target-length", type=float, default=15.0, help="cm of LH2")
    ap.add_argument("--window-ns", type=float, default=100.0, help="GEM integration window for occupancy")
    ap.add_argument("--kmin", type=float, default=0.1, help="minimum photon energy considered [MeV]")
    ap.add_argument("--flat-albedo-kHz-mm2", type=float, default=0.0,
                    help="optional uniform beamline-albedo term added to every plane (unknown; needs Geant4)")
    ap.add_argument("--hadron-charged", type=float, default=1.5, help="charged secondaries per hadronic event")
    ap.add_argument("--hadron-photons", type=float, default=1.0, help="photons (pi0) per hadronic event")
    ap.add_argument("--config-out", default="config/gem_background_rates_analytic.txt")
    ap.add_argument("--no-config", action="store_true")
    args = ap.parse_args()

    E0 = args.beam_energy
    Lt = args.target_length
    n_areal = N_E_LH2 * Lt                         # electrons (protons) per cm^2
    electrons_per_s = args.luminosity / n_areal    # beam electrons per second
    current_uA = electrons_per_s * E_CHARGE * 1e6
    t_target = Lt / LH2_X0
    B_t = bz_axis(TARGET_Z)

    print(f"SoLID DDVCS background estimate  L={args.luminosity:.3g} cm^-2 s^-1")
    print(f"  LH2 target {Lt} cm  -> areal density {n_areal:.3g} e/cm^2, beam current {current_uA:.1f} uA "
          f"({electrons_per_s:.3g} e-/s), t = {t_target:.4f} X0")
    print(f"  Bz(target) = {B_t:.3f} T, Bz(max) = {B_CENTER} T -> magnetic mirror for Moller angles > "
          f"{math.degrees(math.asin(math.sqrt(B_t / B_CENTER))):.1f} deg")

    # total Moller and brems photon fluxes (sanity numbers)
    K = E0 - M_E
    T_grid = np.logspace(math.log10(1.0), math.log10(K / 2), 400)
    sig_moller_1MeV = _trapz([moller_dsigma_dT(T, K) for T in T_grid], T_grid)
    n_moller = electrons_per_s * n_areal * sig_moller_1MeV
    k_grid = np.logspace(math.log10(args.kmin), math.log10(E0 * 0.999), 600)
    n_brems = electrons_per_s * _trapz([brems_spectrum(k, E0, t_target) for k in k_grid], k_grid)
    p_brems = _trapz([brems_spectrum(k, E0, t_target) * k for k in k_grid], k_grid) * electrons_per_s * 1.602e-13
    print(f"  Moller electrons with T > 1 MeV: {n_moller:.3g} /s   (they are confined/mirrored, see below)")
    print(f"  beam bremsstrahlung photons k > {args.kmin} MeV: {n_brems:.3g} /s, power {p_brems:.0f} W, "
          f"all within ~{M_E / E0 * 1e3:.2f} mrad of the beam")

    # hadronic event rate
    sig_eff_virtual = _trapz([equivalent_photon_flux(k, E0) * sigma_gamma_p(k * 1e-3) for k in k_grid], k_grid)
    sig_eff_real = _trapz([brems_spectrum(k, E0, t_target) * N_E_LH2 * (Lt / 2) * sigma_gamma_p(k * 1e-3)
                             for k in k_grid], k_grid) / n_areal
    n_had = electrons_per_s * n_areal * (sig_eff_virtual + sig_eff_real)
    print(f"  hadronic photo/electro-production: sigma_eff = {(sig_eff_virtual + sig_eff_real) / 1e-30:.1f} ub/proton "
          f"-> {n_had:.3g} events/s")

    # Compton-scattered beam-brems photons: dN/dOmega(theta) with scattered energy spectrum
    def compton_dN_dOmega(theta):
        vals, kps = [], []
        for k in k_grid:
            dsig, kp = kn_dsigma_dOmega(k, theta)
            vals.append(brems_spectrum(k, E0, t_target) * N_E_LH2 * (Lt / 2) * dsig)
            kps.append(kp)
        vals = np.array(vals)
        total = _trapz(vals, k_grid) * electrons_per_s
        # interaction-weighted (photon energy dependent) rate
        weighted = _trapz(vals * np.array([gem_interaction_prob(kp) for kp in kps]), k_grid) * electrons_per_s
        return total, weighted

    # ---------------- per plane, per radius
    print()
    header = (f"{'plane':>5} {'R[cm]':>6} {'theta':>6} | {'Moller e':>9} {'Moller-brem g':>13} {'Compton g':>10} "
              f"{'hadron':>8} | {'TOTAL':>9} {'placeholder':>11}   [kHz/mm^2 hits]")
    print(header)
    print("-" * len(header))
    fit_rows = []
    per_plane_totals = {}
    for (pid, zp, rin, rout) in PLANES:
        d = zp - TARGET_Z
        R_grid = np.linspace(rin, rout, 12)
        tot_grid = []
        for R in R_grid:
            theta = math.atan2(R, d)
            proj = math.cos(theta)**3 / d**2        # dOmega -> dA on the plane [1/cm^2]
            # 1. Moller electrons: direct hit only if the helix actually reaches (rin,rout) at this z
            dNdO, T = moller_dN_dOmega(theta, E0, electrons_per_s * n_areal)
            r_reach = helix_reaches(theta, T, zp, rin, rout)
            moller_direct = dNdO * proj if (r_reach is not None and abs(r_reach - R) < 0.05 * R) else 0.0
            # 2. photons radiated by the Moller electron on its way out of the target
            sinth = math.sin(theta)
            L_lh2 = min(TARGET_RADIUS / sinth, Lt / 2 / math.cos(theta))
            t_path = L_lh2 / LH2_X0 + WALL_AL / (AL_X0 * sinth)
            n_gamma = brems_yield(T, t_path, args.kmin)
            # photon energies ~ T/2 typical -> interaction probability at that energy
            p_int = gem_interaction_prob(max(T / 3, args.kmin)) * GEM_SIGNAL_PROB
            moller_brem = dNdO * n_gamma * proj * p_int
            # 3. Compton of beam bremsstrahlung
            c_tot, c_w = compton_dN_dOmega(theta)
            compton = c_w * proj * GEM_SIGNAL_PROB
            # 4. hadronic secondaries (charged: unconfined for p > ~0.1 GeV; photons: pair conversion ~1%)
            had_ang = hadron_angular(theta)
            hadron = n_had * proj * had_ang * (args.hadron_charged * 0.8 + args.hadron_photons * 0.01)
            total = moller_direct + moller_brem + compton + hadron + args.flat_albedo_kHz_mm2 * 1e5
            tot_grid.append(total)
            placeholder = PLACEHOLDER[pid] * 1e5 * (40.0 / R)**2 * args.luminosity / 1.2e37
            if R in (R_grid[0], R_grid[len(R_grid) // 2], R_grid[-1]):
                print(f"{pid:>5} {R:6.1f} {math.degrees(theta):6.1f} | {moller_direct / 1e5:9.4f} {moller_brem / 1e5:13.4f} "
                      f"{compton / 1e5:10.4f} {hadron / 1e5:8.4f} | {total / 1e5:9.4f} {placeholder / 1e5:11.3f}")
        tot_grid = np.array(tot_grid)
        # integrated rate over the plane
        integ = _trapz(tot_grid * 2 * math.pi * R_grid, R_grid)
        integ_ph = _trapz(PLACEHOLDER[pid] * 1e5 * (40.0 / R_grid)**2 * args.luminosity / 1.2e37 * 2 * math.pi * R_grid, R_grid)
        per_plane_totals[pid] = (integ, integ_ph)
        # fit rho0 (40/R)^alpha in log space, area weighted
        w = R_grid
        A = np.vstack([np.ones_like(R_grid), np.log(40.0 / R_grid)]).T
        coef, *_ = np.linalg.lstsq(A * w[:, None], np.log(np.maximum(tot_grid, 1e-30)) * w, rcond=None)
        rho0 = math.exp(coef[0]) / 1e5 * 1.2e37 / args.luminosity   # kHz/mm^2 at reference luminosity
        fit_rows.append((pid, rho0, coef[1]))
        # occupancy at Rin for a strip (0.4 mm x min(rout-rin,100) cm) and a 1 mm^2 pixel
        n_mm2 = tot_grid[0] / 100.0 * args.window_ns * 1e-9   # hits per mm^2 in the window
        occ_strip = n_mm2 * 0.4 * min(rout - rin, 100.0) * 10
        print(f"      plane {pid}: integrated {integ / 1e6:8.1f} MHz (placeholder {integ_ph / 1e6:8.1f} MHz);  "
              f"occupancy@Rin in {args.window_ns:.0f} ns: strip {occ_strip:.3f}, 1 mm^2 pixel {n_mm2:.2e}")

    print()
    print("Fitted rho0 (40/R)^alpha at L = 1.2e37 (for gem:background_file):")
    for pid, rho0, alpha in fit_rows:
        print(f"  plane {pid}: rho0 = {rho0:.4f} kHz/mm^2, alpha = {alpha:.2f}   (placeholder {PLACEHOLDER[pid]}, 2.0)")
    if not args.no_config:
        with open(args.config_out, "w") as f:
            f.write("# Analytic direct-source estimate of the GEM background (scripts/estimate_background_rates.py)\n")
            f.write(f"# L_ref = 1.2e37 cm^-2 s^-1, {Lt} cm LH2, {current_uA:.1f} uA, kmin = {args.kmin} MeV, "
                    f"flat albedo {args.flat_albedo_kHz_mm2} kHz/mm^2\n")
            f.write("# Direct sources only (Moller-brems photons, Compton of beam brems, hadrons); the beamline\n")
            f.write("# albedo term requires a Geant4 beam-on-target simulation and is NOT included.\n")
            f.write("# plane  rho0_kHz_per_mm2  alpha\n")
            for pid, rho0, alpha in fit_rows:
                f.write(f"{pid}  {rho0:.5f}  {alpha:.3f}\n")
        print(f"wrote {args.config_out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
