# GNN tracking for SoLID DDVCS

`tracking:engine=gnn` (built with `-DUSE_ONNX=ON`) finds tracks with a graph neural network that
classifies hit-pair segments, in the style of Exa.TrkX and the ACTS GNN plugin. The model is
trained in PyTorch and runs in C++ through ONNX Runtime, the runtime EICrecon uses. The code is
in `src/gnn/` and `scripts/gnn/`, and the trained model is `config/gnn/solid_gem_in.onnx`.

## Chain

1. **Graph** (`GNNGraph.cc`). Nodes are GEM hits. Edges go from a hit on plane i to a hit on
   plane i+1 or i+2 (the second allows one missed plane), outward in r, and must pass:
   * |Δφ| < 0.25 rad,
   * |z0 − z_target| < 50 cm, where z0 is the straight r–z line through the two hits extrapolated
     to r = 0,
   * |Δφ/Δr| < 1.2 rad/m.
   These cuts keep 99.9% of the e′/μ± segments and about 80% of those of the soft, curling recoil
   proton, which the analysis does not need. The C++ code computes the node features
   (x, y, r, z) and edge features (Δr, Δφ, Δz, z0 offset, Δφ/Δr, plane gap) both for training and
   at reconstruction time.
2. **Edge classifier**: an Interaction Network with node and edge encoders, 4 residual
   message-passing iterations (hidden size 64) and an edge MLP, 160k parameters in total.
   Interface:
   `node_features [N,4], edge_index [2,E], edge_features [E,6] → edge_score [E]`.
3. **Track building** (`gnn:track_builder`):
   * `walkthrough` (default): every path along edges with score > `gnn:edge_threshold` (0.02),
     starting from hits with no incoming edge. Paths are ranked by length and summed score, then
     fitted and accepted greedily without shared hits.
   * `mutual_best`: chains of edges that are the best outgoing edge of their inner hit and the
     best incoming edge of their outer hit. It is faster, but at 10³⁹ a background hit a few
     millimetres from a true hit can win one comparison and cut the track (75% efficiency instead
     of 96%).
4. **Fit**: the native global χ² fitter (RK4 in the field, multiple-scattering covariance,
   beam-line constraint). The hit with the largest pull is removed while χ²/ndf > 10 and more
   than 4 hits remain.
5. **Output**: the same `Track` object as the other engines (truth matching, calorimeter and
   muon-detector flags). An optional calorimeter match is available (`gnn:ec_match_cm`, off by
   default).

## Training

```bash
# 1. export graphs (no tracking is run): one sample per condition, training and validation seeds
jana -Pplugins=solid_ddvcs -Pfastsim:nevents=4000 -Pfastsim:seed=11 -Psolid:luminosity=1e39 \
     -Pgem:readout=pixel -Pgnn:dump=train_pixel_1e39
#    ... likewise for the other conditions and a validation sample with another seed
# 2. train and export to ONNX (PyTorch, onnx, onnxruntime, pandas, numpy)
python3 scripts/gnn/train_gnn.py --train train_pixel_1.2e37 train_pixel_1e38 train_pixel_1e39 train_strip_1.2e37 \
        --val val_pixel_1.2e37 val_pixel_1e38 val_pixel_1e39 val_strip_1.2e37 \
        --out config/gnn/solid_gem_in.onnx --epochs 25
```

The committed model was trained on 16k fast-simulation events (4000 each of pixel readout at
1.2×10³⁷, 10³⁸ and 10³⁹, and strip readout at 1.2×10³⁷): 25 epochs, about 45 minutes on 4 CPU
cores. Validation edge efficiency and purity at score 0.5:

| sample | edge efficiency | edge purity |
|---|---|---|
| pixel 1.2×10³⁷ | 0.9995 | 0.9998 |
| pixel 10³⁸ | 0.998 | 0.994 |
| pixel 10³⁹ | 0.967 | 0.817 |
| strip 1.2×10³⁷ | 0.995 | 0.974 |

ONNX Runtime agrees with PyTorch to 1e-7. `config/gnn/solid_gem_in.json` records the model
configuration and the validation results. The model learns the background of the fast
simulation's placeholder table; retrain it with Geant4 events and a realistic background before
drawing conclusions on real-detector performance.

## Native vs ACTS vs GNN on the same events

These are the same events as in [`acts_tracking.md`](acts_tracking.md); they are independent of
the training samples. Times are wall clock per event with 4 threads and include the fast simulation.

| L [cm⁻²s⁻¹] | readout | engine | events | ε(reco \| geom. acc) | ε(reco \| ≥4 GEM hits) | σ(MM²) core [MeV²] | σ(p)/p μ | σ(θ) μ [mrad] | σ(z_v) [cm] | tracks w/ bkg hit | fake tracks | median χ²/ndf | time / event [ms, 4 threads] |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1.2e+37 | pixel | native | 4000 | 0.825 | 0.998 | 58 | 0.76% | 0.45 | 0.24 | 0.000 | 0.0000 | 0.90 | 8.8 |
| 1.2e+37 | pixel | acts | 4000 | 0.822 | 0.994 | 58 | 0.77% | 0.45 | 0.24 | 0.000 | 0.0000 | 0.78 | 2.6 |
| 1.2e+37 | pixel | gnn | 4000 | 0.822 | 0.994 | 58 | 0.78% | 0.44 | 0.24 | 0.000 | 0.0000 | 0.85 | 2.2 |
| 1.0e+38 | pixel | native | 4000 | 0.784 | 0.981 | 56 | 0.76% | 0.44 | 0.25 | 0.001 | 0.0000 | 0.89 | 10.8 |
| 1.0e+38 | pixel | acts | 4000 | 0.787 | 0.985 | 56 | 0.77% | 0.45 | 0.25 | 0.001 | 0.0000 | 0.78 | 2.6 |
| 1.0e+38 | pixel | gnn | 4000 | 0.789 | 0.987 | 56 | 0.77% | 0.44 | 0.25 | 0.002 | 0.0000 | 0.85 | 2.3 |
| 1.0e+39 | pixel | native | 4000 | 0.766 | 0.940 | 59 | 0.77% | 0.44 | 0.24 | 0.023 | 0.0051 | 0.89 | 80.8 |
| 1.0e+39 | pixel | acts | 4000 | 0.793 | 0.975 | 58 | 0.78% | 0.44 | 0.24 | 0.047 | 0.0342 | 0.80 | 13.0 |
| 1.0e+39 | pixel | gnn | 4000 | 0.780 | 0.959 | 59 | 0.79% | 0.44 | 0.24 | 0.038 | 0.0020 | 0.87 | 7.4 |
| 1.2e+37 | strip | native | 4000 | 0.614 | 0.908 | 52 | 0.76% | 0.45 | 0.25 | 0.006 | 0.0001 | 0.89 | 22.5 |
| 1.2e+37 | strip | acts | 4000 | 0.668 | 0.988 | 53 | 0.77% | 0.45 | 0.25 | 0.006 | 0.0001 | 0.77 | 3.1 |
| 1.2e+37 | strip | gnn | 4000 | 0.659 | 0.975 | 52 | 0.77% | 0.45 | 0.25 | 0.006 | 0.0000 | 0.84 | 2.6 |
| 1.0e+38 | strip | native | 1000 | 0.043 | 0.312 | — | 0.76% | 0.70 | 0.28 | 0.549 | 0.4414 | 2.46 | 525.0 |
| 1.0e+38 | strip | acts | 1000 | 0.029 | 0.250 | — | 0.77% | 0.55 | 0.25 | 0.895 | 0.8828 | 2.24 | 208.9 |
| 1.0e+38 | strip | gnn | 1000 | 0.079 | 0.625 | 52 | 0.84% | 0.62 | 0.30 | 0.652 | 0.5095 | 4.60 | 584.3 |
| 1.0e+39 | strip | native | 200 | 0.000 | 0.000 | — | — | — | — | 1.000 | 1.0000 | 5.51 | 557.1 |
| 1.0e+39 | strip | acts | 200 | 0.000 | 0.000 | — | — | — | — | 1.000 | 1.0000 | 2.74 | 1132.6 |
| 1.0e+39 | strip | gnn | 200 | 0.000 | 0.000 | — | — | — | — | 1.000 | 1.0000 | 7.70 | 1017.8 |

* **Resolution is the same for all three engines**, because the GNN candidates are fitted by the
  native fitter.
* **Efficiency**: the GNN equals the other engines at 1.2×10³⁷ and 10³⁸. At 10³⁹ (pixel) it
  reaches 95.9%, against 94.0% for native and 97.5% for ACTS, while keeping **far fewer fake
  tracks: 0.20%**, against 0.51% for native and 3.4% for ACTS.
* **Saturated strip readout** at 10³⁸, a condition it was not trained on: the GNN reconstructs
  63% of events with three leptons of at least 4 GEM hits, against 31% (native) and 25% (ACTS).
  That is still far from usable: only 8% of geometrically accepted events, with 51% fakes.
* **Speed**: the GNN is the fastest engine below 10³⁹ (2.2–2.6 ms per event, including the fast
  simulation), with one inference per event on a single CPU thread.
