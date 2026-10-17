#!/usr/bin/env python3
"""Train the GNN edge classifier for SoLID GEM track finding and export it to ONNX.

Input: hit graphs written by the JANA chain with -Pgnn:dump=<prefix> (GNNGraphDump_processor):
<prefix>_nodes.csv and <prefix>_edges.csv, with node/edge features and edge truth labels computed
by the same C++ code (src/gnn/GNNGraph.cc) that builds the graphs at reconstruction time.

Model: Interaction Network (Battaglia et al. 2016, as used for HEP tracking by Exa.TrkX):
node and edge encoders, K residual message-passing iterations (edge update from the two end
nodes and the edge; node update from the node and the summed messages on incoming and outgoing
edges), and an edge classifier. The exported ONNX model has
    inputs : node_features float32 [N, 4], edge_index int64 [2, E], edge_features float32 [E, 6]
    output : edge_score    float32 [E]   (sigmoid probability that the edge is a true segment)
and runs in GNNTrackFinder_factory through ONNX Runtime.

usage:
  train_gnn.py --train train_pixel_1e38 train_strip_1.2e37 ... --val val_pixel_1e38 ... \\
               --out config/gnn/solid_gem_in.onnx [--epochs 30 --hidden 64 --iterations 4]
"""
import argparse
import json
import math
import time

import numpy as np
import pandas as pd
import torch
import torch.nn as nn

N_NODE = 4
N_EDGE = 6


def load_graphs(prefix):
    nodes = pd.read_csv(prefix + "_nodes.csv")
    edges = pd.read_csv(prefix + "_edges.csv")
    nf = [f"f{i}" for i in range(N_NODE)]
    ef = [f"f{i}" for i in range(N_EDGE)]
    graphs = []
    eg = dict(tuple(edges.groupby("event")))
    for ev, n in nodes.groupby("event"):
        n = n.sort_values("node")
        e = eg.get(ev)
        if e is None or len(e) == 0:
            continue
        graphs.append(dict(
            x=torch.tensor(n[nf].values, dtype=torch.float32),
            ei=torch.tensor(np.stack([e.src.values, e.dst.values]), dtype=torch.int64),
            ea=torch.tensor(e[ef].values, dtype=torch.float32),
            y=torch.tensor(e.label.values, dtype=torch.float32),
            tag=prefix,
        ))
    return graphs


def mlp(nin, nhid, nout, layers=2, final_act=True):
    mods, d = [], nin
    for _ in range(layers - 1):
        mods += [nn.Linear(d, nhid), nn.LayerNorm(nhid), nn.ReLU()]
        d = nhid
    mods += [nn.Linear(d, nout)]
    if final_act:
        mods += [nn.LayerNorm(nout), nn.ReLU()]
    return nn.Sequential(*mods)


class InteractionNetwork(nn.Module):
    def __init__(self, hidden=64, iterations=4):
        super().__init__()
        self.iterations = iterations
        self.node_enc = mlp(N_NODE, hidden, hidden)
        self.edge_enc = mlp(N_EDGE, hidden, hidden)
        self.edge_net = nn.ModuleList([mlp(3 * hidden, hidden, hidden) for _ in range(iterations)])
        self.node_net = nn.ModuleList([mlp(3 * hidden, hidden, hidden) for _ in range(iterations)])
        self.classifier = mlp(3 * hidden, hidden, 1, layers=3, final_act=False)

    def forward(self, x, edge_index, edge_attr):
        src, dst = edge_index[0], edge_index[1]
        h = self.node_enc(x)
        e = self.edge_enc(edge_attr)
        hidden = h.shape[1]
        idx_dst = dst.unsqueeze(1).expand(-1, hidden)
        idx_src = src.unsqueeze(1).expand(-1, hidden)
        for k in range(self.iterations):
            e = e + self.edge_net[k](torch.cat([h[src], h[dst], e], dim=1))
            agg_in = torch.zeros_like(h).scatter_add(0, idx_dst, e)
            agg_out = torch.zeros_like(h).scatter_add(0, idx_src, e)
            h = h + self.node_net[k](torch.cat([h, agg_in, agg_out], dim=1))
        return self.classifier(torch.cat([h[src], h[dst], e], dim=1)).squeeze(1)


class Exported(nn.Module):
    """ONNX wrapper: probabilities instead of logits."""

    def __init__(self, model):
        super().__init__()
        self.model = model

    def forward(self, node_features, edge_index, edge_features):
        return torch.sigmoid(self.model(node_features, edge_index, edge_features))


def batches(graphs, size, shuffle):
    order = np.random.permutation(len(graphs)) if shuffle else np.arange(len(graphs))
    for i in range(0, len(order), size):
        xs, eis, eas, ys, off = [], [], [], [], 0
        for j in order[i:i + size]:
            g = graphs[j]
            xs.append(g["x"]); eis.append(g["ei"] + off); eas.append(g["ea"]); ys.append(g["y"])
            off += g["x"].shape[0]
        yield torch.cat(xs), torch.cat(eis, dim=1), torch.cat(eas), torch.cat(ys)


@torch.no_grad()
def evaluate(model, graphs, thr=0.5):
    model.eval()
    res = {}
    for tag in sorted({g["tag"] for g in graphs}):
        tp = fp = fn = 0
        for x, ei, ea, y in batches([g for g in graphs if g["tag"] == tag], 64, False):
            p = torch.sigmoid(model(x, ei, ea)) > thr
            t = y > 0.5
            tp += int((p & t).sum()); fp += int((p & ~t).sum()); fn += int((~p & t).sum())
        res[tag] = dict(efficiency=tp / max(1, tp + fn), purity=tp / max(1, tp + fp))
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--train", nargs="+", required=True)
    ap.add_argument("--val", nargs="+", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--epochs", type=int, default=30)
    ap.add_argument("--hidden", type=int, default=64)
    ap.add_argument("--iterations", type=int, default=4)
    ap.add_argument("--batch", type=int, default=32)
    ap.add_argument("--lr", type=float, default=2e-3)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()
    torch.manual_seed(args.seed)
    np.random.seed(args.seed)
    torch.set_num_threads(args.threads)

    t0 = time.time()
    train = [g for p in args.train for g in load_graphs(p)]
    val = [g for p in args.val for g in load_graphs(p)]
    npos = sum(float(g["y"].sum()) for g in train)
    nall = sum(len(g["y"]) for g in train)
    pos_weight = torch.tensor(math.sqrt((nall - npos) / max(npos, 1.0)))
    print(f"loaded {len(train)} train / {len(val)} val graphs in {time.time() - t0:.0f} s; "
          f"{nall} train edges, true fraction {npos / nall:.4f}, pos_weight {float(pos_weight):.2f}", flush=True)

    model = InteractionNetwork(args.hidden, args.iterations)
    print(f"Interaction Network: hidden {args.hidden}, {args.iterations} iterations, "
          f"{sum(p.numel() for p in model.parameters())} parameters", flush=True)
    opt = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=1e-4)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=args.epochs)
    lossf = nn.BCEWithLogitsLoss(pos_weight=pos_weight)
    best, best_state = -1.0, None
    for epoch in range(args.epochs):
        model.train()
        tl, nb, t1 = 0.0, 0, time.time()
        for x, ei, ea, y in batches(train, args.batch, True):
            opt.zero_grad()
            loss = lossf(model(x, ei, ea), y)
            loss.backward()
            opt.step()
            tl += float(loss); nb += 1
        sched.step()
        res = evaluate(model, val)
        # model selection: worst-case F1 over the validation configurations
        score = min(2 * r["efficiency"] * r["purity"] / max(1e-9, r["efficiency"] + r["purity"]) for r in res.values())
        line = " ".join(f"{k.split('/')[-1]}: eff {r['efficiency']:.4f} pur {r['purity']:.4f}" for k, r in res.items())
        print(f"epoch {epoch:2d} loss {tl / nb:.4f} ({time.time() - t1:.0f} s)  {line}", flush=True)
        if score > best:
            best, best_state = score, {k: v.clone() for k, v in model.state_dict().items()}
    model.load_state_dict(best_state)
    model.eval()
    final = evaluate(model, val)

    # export to ONNX and check it against PyTorch
    g = val[0]
    exp = Exported(model).eval()
    torch.onnx.export(exp, (g["x"], g["ei"], g["ea"]), args.out, dynamo=False, opset_version=17,
                      input_names=["node_features", "edge_index", "edge_features"], output_names=["edge_score"],
                      dynamic_axes={"node_features": {0: "nodes"}, "edge_index": {1: "edges"},
                                    "edge_features": {0: "edges"}, "edge_score": {0: "edges"}})
    import onnxruntime as ort
    sess = ort.InferenceSession(args.out, providers=["CPUExecutionProvider"])
    maxdiff = 0.0
    for g in val[:50]:
        ref = exp(g["x"], g["ei"], g["ea"]).detach().numpy()
        out = sess.run(None, {"node_features": g["x"].numpy(), "edge_index": g["ei"].numpy(),
                              "edge_features": g["ea"].numpy()})[0]
        maxdiff = max(maxdiff, float(np.abs(out - ref).max()))
    print(f"exported {args.out}; max |ONNX Runtime - PyTorch| over 50 graphs = {maxdiff:.2e}")
    meta = dict(model="InteractionNetwork", hidden=args.hidden, iterations=args.iterations,
                node_features=["x/100", "y/100", "r/100", "z/100"],
                edge_features=["dr/100", "dphi", "dz/100", "(z0-z_target)/100", "dphi/dr*100 [rad/m]", "plane gap - 1"],
                train=args.train, val=args.val, epochs=args.epochs, edge_threshold=0.5,
                validation=final, onnx_vs_torch_maxdiff=maxdiff)
    with open(args.out.replace(".onnx", ".json"), "w") as f:
        json.dump(meta, f, indent=2)
    print(json.dumps(final, indent=1))


if __name__ == "__main__":
    main()
