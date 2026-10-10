#!/usr/bin/env python3
"""Write a HepMC3 ASCII file laid out exactly as EpIC writes DDVCS events, for testing HepMC3Reader.

The event construction follows EpIC src/modules/kinematic/DDVCS/DDVCSKinematicDefault.cpp step by
step (lepton beam along -z, target with E = 1.00001 Mp along +z, ROOT rotation/boost conventions,
lepton pair built from thetaL/phiL around the recoil proton in the dilepton rest frame), and the
file layout follows HepMC3::WriterAscii as used by EpIC WriterHepMC3 (status 4 beams, 1 final,
13 incoming gamma*, 3 outgoing gamma*; single-parent vertices written implicitly; run information
rewritten after the last event). The generated (Q2, y, t, Q'2, phi, thetaL, phiL, phiS) are written
to a CSV next to it so the test can compare the kinematics recomputed by the reader.

usage: make_epic_like_hepmc.py out.hepmc out_truth.csv [nevents] [seed]
"""
import math
import random
import sys

ME, MMU, MP = 0.000510999, 0.1056583755, 0.938272088


# --- minimal TVector3 / TLorentzVector with ROOT conventions -------------------------------------
def add(a, b): return [x + y for x, y in zip(a, b)]
def sub(a, b): return [x - y for x, y in zip(a, b)]
def mul(a, s): return [x * s for x in a]
def dot(a, b): return sum(x * y for x, y in zip(a, b))
def cross(a, b): return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
def mag(a): return math.sqrt(dot(a, a))
def unit(a): m = mag(a); return mul(a, 1 / m) if m > 0 else a


def rotate(v, angle, axis):
    """TVector3::Rotate(angle, axis): right-handed rotation (Rodrigues)."""
    k = unit(axis)
    c, s = math.cos(angle), math.sin(angle)
    return add(add(mul(v, c), mul(cross(k, v), s)), mul(k, dot(k, v) * (1 - c)))


def rotY(v, a): return rotate(v, a, [0, 1, 0])
def rotZ(v, a): return rotate(v, a, [0, 0, 1])


class LV:
    def __init__(self, p, E): self.p, self.E = list(p), E
    @staticmethod
    def dirE(d, E, m): pm = math.sqrt(max(0.0, E * E - m * m)); return LV(mul(unit(d), pm), E)
    @staticmethod
    def pm(p, m): return LV(p, math.sqrt(dot(p, p) + m * m))
    def __add__(self, o): return LV(add(self.p, o.p), self.E + o.E)
    def __sub__(self, o): return LV(sub(self.p, o.p), self.E - o.E)
    def m2(self): return self.E ** 2 - dot(self.p, self.p)
    def beta(self): return mul(self.p, 1 / self.E)
    def boost(self, b):            # TLorentzVector::Boost(b)
        b2 = dot(b, b)
        if b2 == 0: return LV(self.p, self.E)
        g = 1 / math.sqrt(1 - b2)
        bp = dot(b, self.p)
        g2 = (g - 1) / b2
        return LV(add(self.p, mul(b, g2 * bp + g * self.E)), g * (self.E + bp))
    def rot(self, angle, axis): return LV(rotate(self.p, angle, axis), self.E)


def epic_event(Ee, Ep, y, Q2, t, Q2p, phi, thetaL, phiL, phiS, mlep):
    e_LAB = LV.dirE([0, 0, -1], Ee, ME)
    p_LAB = LV.dirE([0, 0, 1], Ep, MP)
    b_LT = p_LAB.beta()
    e_TAR, p_TAR = e_LAB.boost(mul(b_LT, -1)), p_LAB.boost(mul(b_LT, -1))
    E_eS = e_TAR.E * (1 - y)
    cos_te = 1 - Q2 / (2 * e_TAR.E * E_eS)
    eS_TAR = LV.dirE(e_TAR.p, E_eS, ME)
    eS_TAR = LV(rotY(eS_TAR.p, math.acos(cos_te)), eS_TAR.E)
    g_TAR = e_TAR - eS_TAR
    s4 = (g_TAR + p_TAR).m2()
    if (s4 - Q2p - MP ** 2) ** 2 - 4 * Q2p * MP ** 2 < 0: return None
    m1, m2, m3, m4 = -Q2, MP ** 2, Q2p, MP ** 2
    E1 = (s4 + m1 - m2) / (2 * math.sqrt(s4)); E3 = (s4 + m3 - m4) / (2 * math.sqrt(s4))
    p1 = math.sqrt(E1 ** 2 - m1); p3 = math.sqrt(E3 ** 2 - m3)
    c0 = ((m1 - m3 - m2 + m4) / (2 * math.sqrt(s4))) ** 2
    tmin, tmax = c0 - (p1 - p3) ** 2, c0 - (p1 + p3) ** 2
    if t > tmin or t < tmax: return None
    b_TG = (g_TAR + p_TAR).beta()
    g_GP, p_GP = g_TAR.boost(mul(b_TG, -1)), p_TAR.boost(mul(b_TG, -1))
    pex = math.sqrt(((s4 - Q2p - MP ** 2) ** 2 - 4 * Q2p * MP ** 2) / (4 * s4))
    ex_GP = LV(mul(unit(g_GP.p), pex), math.sqrt(pex ** 2 + Q2p))
    pS_GP = LV.pm(mul(ex_GP.p, -1), MP)
    cth = (t - (p_GP.E - pS_GP.E) ** 2 + dot(p_GP.p, p_GP.p) + dot(pS_GP.p, pS_GP.p)) / (2 * mag(p_GP.p) * mag(pS_GP.p))
    if abs(cth) > 1: return None
    ex_GP = LV(rotY(ex_GP.p, math.acos(cth)), ex_GP.E)
    pS_GP = LV(rotY(pS_GP.p, math.acos(cth)), pS_GP.E)
    ex_TAR, pS_TAR = ex_GP.boost(b_TG), pS_GP.boost(b_TG)
    qax = unit(g_TAR.p)
    pS_TAR, ex_TAR = pS_TAR.rot(phi, qax), ex_TAR.rot(phi, qax)
    b_TE = ex_TAR.beta()
    p_EXC, pS_EXC = p_TAR.boost(mul(b_TE, -1)), pS_TAR.boost(mul(b_TE, -1))
    El = math.sqrt(Q2p) / 2
    l1 = LV.dirE(pS_EXC.p, El, mlep); l2 = LV.dirE(mul(pS_EXC.p, -1), El, mlep)
    ax = unit(cross(p_EXC.p, pS_EXC.p))
    l1, l2 = l1.rot(thetaL, ax), l2.rot(thetaL, ax)
    l1, l2 = l1.rot(phiL, unit(pS_EXC.p)), l2.rot(phiL, unit(pS_EXC.p))
    l1, l2 = l1.boost(b_TE), l2.boost(b_TE)
    out = {}
    for name, v in (("g", g_TAR), ("ex", ex_TAR), ("eS", eS_TAR), ("pS", pS_TAR), ("l1", l1), ("l2", l2)):
        out[name] = LV(rotZ(v.p, phiS), v.E).boost(b_LT)
    out["e"], out["p"] = e_LAB, p_LAB
    return out


def fmt(x): return "%.16e" % x


def main():
    out, truth = sys.argv[1], sys.argv[2]
    n = int(sys.argv[3]) if len(sys.argv) > 3 else 200
    rng = random.Random(int(sys.argv[4]) if len(sys.argv) > 4 else 1)
    Ee, Ep = 11.0, 1.00001 * MP
    with open(out, "w") as f, open(truth, "w") as ft:
        f.write("HepMC::Version 3.02.06\nHepMC::Asciiv3-START\nW weight_default\n")
        f.write("T EpIC|1.1.6|test fixture\n")
        ft.write("event,Q2,y,t,Q2p,phi,thetaL,phiL,phiS\n")
        i = 0
        while i < n:
            Q2 = rng.uniform(1.0, 4.0); y = rng.uniform(0.2, 0.8); Q2p = rng.uniform(2.0, 6.0)
            t = -rng.uniform(0.05, 1.0)
            phi, thetaL = rng.uniform(0.05, 6.23), rng.uniform(0.79, 2.35)
            phiL, phiS = rng.uniform(0.05, 6.23), rng.uniform(0, 2 * math.pi)
            ev = epic_event(Ee, Ep, y, Q2, t, Q2p, phi, thetaL, phiL, phiS, MMU)
            if ev is None: continue
            mass = {"e": ME, "p": MP, "eS": ME, "g": -math.sqrt(Q2), "pS": MP, "ex": math.sqrt(Q2p), "l1": MMU, "l2": MMU}
            # particle ids as WriterAscii assigns them: e(1) p(2) eS(3) g(4) pS(5) ex(6) mu-(7) mu+(8)
            f.write("E %d 2 8\nU GEV MM\nW 1.0000000000000000e+00\n" % i)
            f.write("A 0 GenCrossSection 1.0000000000000000e+00 0.0 -1 -1\n")
            P = lambda pid, par, pdg, k, st: f.write("P %d %d %d %s %s %s %s %s %d\n" % (
                pid, par, pdg, fmt(ev[k].p[0]), fmt(ev[k].p[1]), fmt(ev[k].p[2]), fmt(ev[k].E), fmt(mass[k]), st))
            P(1, 0, 11, "e", 4)
            P(2, 0, 2212, "p", 4)
            P(3, 1, 11, "eS", 1)           # single-parent vertex: parent given as particle id
            P(4, 1, 22, "g", 13)
            f.write("V -2 0 [4,2]\n")
            P(5, -2, 2212, "pS", 1)
            P(6, -2, 22, "ex", 3)
            P(7, 6, 13, "l1", 1)
            P(8, 6, -13, "l2", 1)
            ft.write("%d,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g\n" % (i, Q2, y, t, Q2p, phi, thetaL, phiL, phiS))
            i += 1
        # EpIC rewrites the run information after the last event
        f.write("W weight_default\nT EpIC|1.1.6|test fixture\n")
        f.write("A integrated_cross_section_value 1.2345000000000000e-02\n")
        f.write("A integrated_cross_section_uncertainty 1.0000000000000000e-04\n")
        f.write("HepMC::Asciiv3-END\n")


if __name__ == "__main__":
    main()
