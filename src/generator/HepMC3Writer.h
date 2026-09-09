// Minimal HepMC3 ASCII (Asciiv3) writer for the DDVCS generator output, readable by
// DDG4's HEPMC3FileReader (scripts/run_ddsim_ddvcs.py) and by HepMC3 itself.
#pragma once

#include <fstream>
#include <string>
#include <vector>
#include "datamodel/SoLIDDataModel.h"

namespace solid {

class HepMC3Writer {
public:
  explicit HepMC3Writer(const std::string& filename);
  ~HepMC3Writer();
  bool good() const { return m_out.good(); }
  /// Write one event. The particle record must follow the DDVCSGenerator layout
  /// (0: beam e-, 1: target p, 2: e', 3: p', 4: gamma*, 5: mu-, 6: mu+).
  void write(int event_number, const std::vector<MCParticle>& particles, double weight = 1.0);

private:
  std::ofstream m_out;
};

}  // namespace solid
