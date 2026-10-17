// ONNX Runtime inference of the GNN edge classifier (model trained and exported by
// scripts/gnn/train_gnn.py). Model interface:
//   inputs : node_features float32 [N, 4], edge_index int64 [2, E], edge_features float32 [E, 6]
//   output : edge_score    float32 [E]
#pragma once

#include <memory>
#include <string>
#include <vector>
#include "gnn/GNNGraph.h"

namespace Ort {
struct Env;
struct Session;
}  // namespace Ort

namespace solid {

class GNNEdgeClassifier {
public:
  /// Loads the model; one intra-op thread (JANA runs events in parallel).
  explicit GNNEdgeClassifier(const std::string& model_path);
  ~GNNEdgeClassifier();
  /// Edge scores in [0, 1], one per edge of the graph.
  std::vector<float> score(const GNNGraph& graph) const;

private:
  std::unique_ptr<Ort::Session> m_session;
};

}  // namespace solid
