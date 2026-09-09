// solid_ddvcs_gen: stand-alone DDVCS generator writing HepMC3 ASCII for ddsim/DDG4.
//   solid_ddvcs_gen -n 10000 -o ddvcs.hepmc [-s seed] [-E 11] [--target-z -315]
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include "DDVCSGenerator.h"
#include "HepMC3Writer.h"

int main(int argc, char** argv) {
  int n = 1000;
  unsigned seed = 12345;
  std::string out = "ddvcs.hepmc";
  solid::DDVCSConfig cfg;
  for (int i = 1; i < argc; ++i) {
    auto next = [&](double& v) { if (i + 1 < argc) v = std::atof(argv[++i]); };
    if (!std::strcmp(argv[i], "-n") && i + 1 < argc) n = std::atoi(argv[++i]);
    else if (!std::strcmp(argv[i], "-s") && i + 1 < argc) seed = static_cast<unsigned>(std::atol(argv[++i]));
    else if (!std::strcmp(argv[i], "-o") && i + 1 < argc) out = argv[++i];
    else if (!std::strcmp(argv[i], "-E")) next(cfg.beam_energy);
    else if (!std::strcmp(argv[i], "--target-z")) next(cfg.target_z);
    else if (!std::strcmp(argv[i], "--q2-min")) next(cfg.Q2_min);
    else if (!std::strcmp(argv[i], "--q2-max")) next(cfg.Q2_max);
    else if (!std::strcmp(argv[i], "--qp2-min")) next(cfg.Qp2_min);
    else if (!std::strcmp(argv[i], "--qp2-max")) next(cfg.Qp2_max);
    else if (!std::strcmp(argv[i], "--tabs-max")) next(cfg.tabs_max);
    else {
      std::cerr << "usage: solid_ddvcs_gen -n N -o file.hepmc [-s seed] [-E Ebeam] [--target-z z] "
                   "[--q2-min/--q2-max/--qp2-min/--qp2-max/--tabs-max v]\n";
      return 1;
    }
  }
  solid::DDVCSGenerator gen(cfg, seed);
  solid::HepMC3Writer writer(out);
  if (!writer.good()) { std::cerr << "cannot open " << out << "\n"; return 2; }
  std::vector<solid::MCParticle> parts;
  solid::DDVCSKinematics kin;
  double mm2_max_dev = 0;
  for (int i = 0; i < n; ++i) {
    if (!gen.generate(parts, kin)) { std::cerr << "generation failed\n"; return 3; }
    writer.write(i + 1, parts);
    double mm2 = solid::DDVCSGenerator::missingMass2(cfg.beam_energy, parts[2].p4, parts[6].p4, parts[5].p4);
    mm2_max_dev = std::max(mm2_max_dev, std::fabs(mm2 - solid::kProtonMass * solid::kProtonMass));
  }
  std::cout << "wrote " << n << " DDVCS events to " << out << " (max |MM2 - Mp2| = " << mm2_max_dev << " GeV^2)\n";
  return 0;
}
