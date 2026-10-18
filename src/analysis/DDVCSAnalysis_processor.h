// Event-level output for the DDVCS missing-mass study: a CSV ntuple and a JSON summary with
// acceptance, reconstruction efficiency, track-quality and missing-mass resolution figures.
#pragma once

#include <JANA/JEventProcessor.h>
#include <fstream>
#include <string>
#include <vector>
#include "datamodel/SoLIDDataModel.h"
#include "services/SoLIDGeometry_service.h"

namespace solid {

class DDVCSAnalysis_processor : public JEventProcessor {
public:
  DDVCSAnalysis_processor();
  void Init() override;
  void ProcessSequential(const JEvent& event) override;
  void Finish() override;

private:
  struct Stats {
    std::vector<double> v;
    void add(double x) { v.push_back(x); }
    size_t n() const { return v.size(); }
    double mean() const;
    double rms() const;
    double quantile(double q) const;   // requires sorted copy
    double sigma68() const { return 0.5 * (quantile(0.84135) - quantile(0.15865)); }
    double sigma95() const { return 0.5 * (quantile(0.97725) - quantile(0.02275)); }
    double gaussCoreSigma(double& mean_out) const;  // iterative 2.5-sigma clipped Gaussian estimate
  };

  Input<DDVCSEvent> m_ev_in{this};
  Input<GenKinematics> m_gen_in{this};
  Input<Track> m_tracks_in{this};
  Input<MCParticle> m_mc_in{this};
  Service<SoLIDGeometry_service> m_geo{this};

  Parameter<std::string> m_csv{this, "analysis:output", "ddvcs_events.csv", "Per-event CSV output"};
  Parameter<std::string> m_json{this, "analysis:summary", "ddvcs_summary.json", "JSON summary output"};
  Parameter<std::string> m_label{this, "analysis:label", "", "Free-text label stored in the summary"};

  std::ofstream m_out;
  long m_n_events = 0, m_n_geom = 0, m_n_found_geom = 0, m_n_accepted = 0, m_n_found = 0, m_n_found_accepted = 0, m_n_all_matched = 0;
  long m_n_tracks = 0, m_n_tracks_with_bkg = 0, m_n_fake_tracks = 0;
  // per lepton in geometrically accepted events with >= 4 GEM hits: [e, mu+, mu-]
  long m_lep_total[3] = {0, 0, 0}, m_lep_tracked[3] = {0, 0, 0}, m_lep_pid_ok[3] = {0, 0, 0};
  Stats m_mm2, m_mm, m_mm2_clean, m_dpp_e, m_dpp_mu, m_dth_e, m_dth_mu, m_dph_mu, m_dQ2, m_dQp2, m_dt, m_chi2ndf, m_dzv;
  Stats m_dphi, m_dcos_thl, m_dphi_l;
  double m_lumi = 0;
  double m_weight_sum = 0, m_weight_sum_found = 0;
  int m_gen_source = 0;
};

}  // namespace solid
