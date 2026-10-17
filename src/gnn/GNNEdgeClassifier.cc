#include "GNNEdgeClassifier.h"
#include <onnxruntime_cxx_api.h>
#include <array>
#include <stdexcept>

namespace solid {

namespace {
Ort::Env& ortEnv() {
  static Ort::Env env(ORT_LOGGING_LEVEL_WARNING, "solid_gnn");
  return env;
}
}  // namespace

GNNEdgeClassifier::GNNEdgeClassifier(const std::string& model_path) {
  Ort::SessionOptions opts;
  opts.SetIntraOpNumThreads(1);
  opts.SetInterOpNumThreads(1);
  opts.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_ENABLE_ALL);
  m_session = std::make_unique<Ort::Session>(ortEnv(), model_path.c_str(), opts);
}

GNNEdgeClassifier::~GNNEdgeClassifier() = default;

std::vector<float> GNNEdgeClassifier::score(const GNNGraph& g) const {
  const auto N = static_cast<std::int64_t>(g.numNodes());
  const auto E = static_cast<std::int64_t>(g.numEdges());
  if (E == 0) return {};
  auto mem = Ort::MemoryInfo::CreateCpu(OrtArenaAllocator, OrtMemTypeDefault);
  std::vector<std::int64_t> edge_index(2 * E);
  std::copy(g.src.begin(), g.src.end(), edge_index.begin());
  std::copy(g.dst.begin(), g.dst.end(), edge_index.begin() + E);
  const std::array<std::int64_t, 2> nshape{N, GNNGraph::kNodeFeatures}, ishape{2, E}, eshape{E, GNNGraph::kEdgeFeatures};
  std::array<Ort::Value, 3> inputs{
      Ort::Value::CreateTensor<float>(mem, const_cast<float*>(g.node_features.data()), g.node_features.size(),
                                      nshape.data(), nshape.size()),
      Ort::Value::CreateTensor<std::int64_t>(mem, edge_index.data(), edge_index.size(), ishape.data(), ishape.size()),
      Ort::Value::CreateTensor<float>(mem, const_cast<float*>(g.edge_features.data()), g.edge_features.size(),
                                      eshape.data(), eshape.size())};
  const std::array<const char*, 3> in_names{"node_features", "edge_index", "edge_features"};
  const std::array<const char*, 1> out_names{"edge_score"};
  auto out = m_session->Run(Ort::RunOptions{nullptr}, in_names.data(), inputs.data(), inputs.size(), out_names.data(),
                            out_names.size());
  const float* s = out[0].GetTensorData<float>();
  if (out[0].GetTensorTypeAndShapeInfo().GetElementCount() != static_cast<size_t>(E))
    throw std::runtime_error("GNN model returned a wrong number of edge scores");
  return std::vector<float>(s, s + E);
}

}  // namespace solid
