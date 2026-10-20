// Minimal HepMC3 ASCII reader for external DDVCS generators, in particular EpIC
// (https://github.com/pawelsznajder/epic, DDVCSGeneratorService with WriterHepMC3 "ascii").
//
// Each event is converted to the internal event record used by the fast simulation and written
// by HepMC3Writer:
//    0 beam lepton (status 4)    1 target (status 4)    2 scattered lepton (status 1)
//    3 recoil proton (status 1)  4 outgoing gamma* (status 11)   5 l- (status 1)   6 l+ (status 1)
//    7... any other final-state particle (e.g. radiated photons), status 1
//
// Frame: SoLID has the beam along +z. EpIC uses the collider convention with the lepton beam
// along -z, so when the beam lepton has p_z < 0 the event is rotated by pi around the y axis
// (x -> -x, z -> -z; a proper rotation, so all angles and kinematics are unchanged).
// Units are taken from the "U" line (GEV/MEV, MM/CM) and converted to GeV and cm.
// Vertex: EpIC writes no vertex positions; vertex_in_file() is then false and the caller places
// the event in the target.
#pragma once

#include <fstream>
#include <map>
#include <string>
#include <vector>
#include "datamodel/SoLIDDataModel.h"

namespace solid {

struct HepMC3RunInfo {
  std::string tool;                    // "T" line(s): generator name and version
  std::map<std::string, std::string> attributes;  // run attributes ("A name value")
  double xsec = 0, xsec_err = 0;       // integrated cross section, if provided (generator units)
  bool has_xsec = false;
};

class HepMC3Reader {
public:
  explicit HepMC3Reader(const std::string& filename);
  bool good() const { return m_in.good() || m_have_pending; }

  /// Read the next event. Returns false at the end of the file. Throws std::runtime_error if
  /// the event does not have the DDVCS topology (beam lepton, target, l', p', l- l+ pair).
  bool next(std::vector<MCParticle>& particles);

  int eventNumber() const { return m_event_number; }
  double weight() const { return m_weight; }
  /// Per-event "GenCrossSection" attribute (EpIC stores the generator weight there), 0 if absent.
  double eventCrossSection() const { return m_event_xsec; }
  bool vertexInFile() const { return m_vertex_in_file; }
  bool flipped() const { return m_flipped; }
  /// Run information seen so far. EpIC writes the integrated cross section after the last event,
  /// so it is complete only once next() has returned false.
  const HepMC3RunInfo& runInfo() const { return m_run; }

private:
  struct RawParticle {
    int id = 0, parent = 0, pdg = 0, status = 0;
    double px = 0, py = 0, pz = 0, e = 0, m = 0;
  };
  struct RawVertex {
    std::vector<int> in;
    double x = 0, y = 0, z = 0;
    bool has_pos = false;
  };
  bool readEvent();
  void parseRunLine(const std::string& line);

  std::ifstream m_in;
  std::string m_pending;     // "E" line already read that starts the next event
  bool m_have_pending = false;
  HepMC3RunInfo m_run;

  // current event
  int m_event_number = 0;
  double m_weight = 1.0;
  double m_event_xsec = 0;
  double m_mom_unit = 1.0, m_len_unit = 0.1;   // to GeV, cm
  std::vector<RawParticle> m_parts;
  std::map<int, RawVertex> m_vertices;
  bool m_vertex_in_file = false;
  bool m_flipped = false;
};

}  // namespace solid
