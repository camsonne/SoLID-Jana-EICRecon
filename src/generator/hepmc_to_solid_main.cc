// solid_hepmc_to_solid: convert HepMC3 ASCII events from an external DDVCS generator (EpIC) to the
// SoLID frame for ddsim/DDG4: beam along +z (EpIC has it along -z), vertex placed uniformly along
// the LH2 target, internal event-record layout of solid_ddvcs_gen (HepMC3Reader.h).
//   solid_hepmc_to_solid -i epic_out_file.txt -o ddvcs_solid.hepmc [-n N] [-s seed] [--target-z -315] [--target-half 7.5]
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <random>
#include <stdexcept>
#include "HepMC3Reader.h"
#include "HepMC3Writer.h"

int main(int argc, char** argv) {
  std::string in, out = "ddvcs_solid.hepmc";
  long nmax = 0;
  unsigned seed = 12345;
  double target_z = -315.0, target_half = 7.5;   // cm, as in SoLIDGeometry_service
  bool keep_vertex = false;
  for (int i = 1; i < argc; ++i) {
    auto next = [&](double& v) { if (i + 1 < argc) v = std::atof(argv[++i]); };
    if (!std::strcmp(argv[i], "-i") && i + 1 < argc) in = argv[++i];
    else if (!std::strcmp(argv[i], "-o") && i + 1 < argc) out = argv[++i];
    else if (!std::strcmp(argv[i], "-n") && i + 1 < argc) nmax = std::atol(argv[++i]);
    else if (!std::strcmp(argv[i], "-s") && i + 1 < argc) seed = static_cast<unsigned>(std::atol(argv[++i]));
    else if (!std::strcmp(argv[i], "--target-z")) next(target_z);
    else if (!std::strcmp(argv[i], "--target-half")) next(target_half);
    else if (!std::strcmp(argv[i], "--keep-vertex")) keep_vertex = true;
    else {
      std::cerr << "usage: solid_hepmc_to_solid -i in.hepmc -o out.hepmc [-n N] [-s seed] [--target-z z] "
                   "[--target-half h] [--keep-vertex]\n";
      return 1;
    }
  }
  if (in.empty()) { std::cerr << "no input file (-i)\n"; return 1; }
  solid::HepMC3Reader reader(in);
  if (!reader.good()) { std::cerr << "cannot open " << in << "\n"; return 2; }
  solid::HepMC3Writer writer(out);
  if (!writer.good()) { std::cerr << "cannot open " << out << "\n"; return 2; }

  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> U(-1.0, 1.0);
  std::vector<solid::MCParticle> parts;
  long n = 0;
  try {
    while ((nmax <= 0 || n < nmax) && reader.next(parts)) {
      if (!keep_vertex || !reader.vertexInFile()) {
        const solid::Vec3 vtx(0, 0, target_z + target_half * U(rng));
        for (auto& p : parts) p.vertex = vtx;
      }
      writer.write(static_cast<int>(++n), parts, reader.weight());
    }
  } catch (const std::exception& e) {
    std::cerr << e.what() << "\n";
    return 3;
  }
  const auto& run = reader.runInfo();
  std::cout << "converted " << n << " events from " << in << " to " << out
            << (reader.flipped() ? " (rotated: beam was along -z)" : "") << "\n";
  if (!run.tool.empty()) std::cout << "generator: " << run.tool << "\n";
  if (run.has_xsec) std::cout << "integrated cross section: " << run.xsec << " +- " << run.xsec_err << " (generator units)\n";
  return 0;
}
