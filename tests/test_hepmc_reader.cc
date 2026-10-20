// Reads an EpIC-like HepMC3 file (tests/make_epic_like_hepmc.py) with HepMC3Reader and checks that
// the frame conversion and DDVCSKinematicsCalc reproduce the generated EpIC variables.
//   test_hepmc_reader events.hepmc truth.csv
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "generator/DDVCSKinematicsCalc.h"
#include "generator/HepMC3Reader.h"

using namespace solid;

namespace {
int g_fail = 0;
void check(bool ok, const std::string& what) {
  if (!ok) { ++g_fail; if (g_fail < 20) std::cerr << "FAIL: " << what << "\n"; }
}
double angDiff(double a, double b) { return std::fabs(std::remainder(a - b, 2 * M_PI)); }
}  // namespace

int main(int argc, char** argv) {
  if (argc < 3) { std::cerr << "usage: test_hepmc_reader events.hepmc truth.csv\n"; return 2; }
  HepMC3Reader reader(argv[1]);
  std::ifstream truth(argv[2]);
  std::string line;
  std::getline(truth, line);  // header
  std::vector<MCParticle> parts;
  int n = 0;
  double max_dphi = 0, max_dphil = 0, max_dcos = 0, max_dQ2 = 0, max_dt = 0, max_dQp2 = 0, max_dy = 0, max_dmm2 = 0;
  while (reader.next(parts)) {
    if (!std::getline(truth, line)) { check(false, "more events than truth rows"); break; }
    std::istringstream is(line);
    std::vector<double> v;
    std::string item;
    while (std::getline(is, item, ',')) v.push_back(std::stod(item));
    const double Q2 = v[1], y = v[2], t = v[3], Qp2 = v[4], phi = v[5], thetaL = v[6], phiL = v[7];
    check(parts.size() == 7, "7 particles");
    check(parts[0].p4.p.z > 0, "beam along +z after conversion");
    check(parts[0].pdg == 11 && parts[2].pdg == 11 && parts[3].pdg == 2212, "e, e', p' identified");
    check(parts[5].pdg == 13 && parts[6].pdg == -13, "mu- mu+ identified");
    check(parts[4].status == 11 && parts[5].parent == 4, "gamma* record");
    check(std::fabs(parts[5].charge + 1) < 1e-9 && std::fabs(parts[6].charge - 1) < 1e-9, "muon charges");
    const auto k = computeDDVCSKinematics(parts[0].p4, parts[1].p4, parts[2].p4, parts[5].p4, parts[6].p4);
    max_dQ2 = std::max(max_dQ2, std::fabs(k.Q2 - Q2));
    max_dy = std::max(max_dy, std::fabs(k.y - y));
    max_dt = std::max(max_dt, std::fabs(k.t - t));
    max_dQp2 = std::max(max_dQp2, std::fabs(k.Qp2 - Qp2));
    max_dmm2 = std::max(max_dmm2, std::fabs(k.mm2 - kProtonMass * kProtonMass));
    max_dcos = std::max(max_dcos, std::fabs(k.cos_theta_l - std::cos(thetaL)));
    max_dphil = std::max(max_dphil, angDiff(k.phi_l, phiL));
    max_dphi = std::max(max_dphi, angDiff(k.phi_trento, phi));
    ++n;
  }
  std::printf("events %d  max|dQ2| %.2e  max|dy| %.2e  max|dt| %.2e  max|dQ'2| %.2e  max|dMM2| %.2e\n", n, max_dQ2,
              max_dy, max_dt, max_dQp2, max_dmm2);
  std::printf("max|dcos(thetaL)| %.2e  max|dphiL| %.2e  max|dphi| %.2e\n", max_dcos, max_dphil, max_dphi);
  check(n > 0, "events read");
  check(reader.flipped(), "EpIC frame detected and rotated");
  check(!reader.vertexInFile(), "EpIC writes no vertex");
  // EpIC neglects the electron mass in the scattering angle: |dQ2| ~ 1e-6 GeV^2
  check(max_dQ2 < 1e-5 && max_dy < 1e-6 && max_dt < 1e-6 && max_dQp2 < 1e-6 && max_dmm2 < 1e-6, "invariants");
  check(max_dcos < 1e-6 && max_dphil < 1e-6, "lepton-pair angles match EpIC thetaL, phiL");
  check(max_dphi < 1e-6, "Trento phi matches EpIC phi");
  const auto& run = reader.runInfo();
  check(run.has_xsec && std::fabs(run.xsec - 1.2345e-2) < 1e-12 && std::fabs(run.xsec_err - 1e-4) < 1e-12,
        "integrated cross section read from the trailing run information");
  check(std::fabs(reader.weight() - 1.0) < 1e-12, "event weight not overwritten by the trailing 'W weight_default'");
  check(run.tool.find("EpIC") != std::string::npos, "generator name");
  if (g_fail) { std::cerr << g_fail << " check(s) failed\n"; return 1; }
  std::cout << "OK\n";
  return 0;
}
