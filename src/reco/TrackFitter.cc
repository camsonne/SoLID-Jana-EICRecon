#include "TrackFitter.h"
#include <algorithm>
#include <cmath>

namespace solid {

bool invertMatrix(std::vector<double>& A, int n) {
  std::vector<double> I(n * n, 0.0);
  for (int i = 0; i < n; ++i) I[i * n + i] = 1.0;
  for (int c = 0; c < n; ++c) {
    int piv = c;
    double best = std::fabs(A[c * n + c]);
    for (int r = c + 1; r < n; ++r) {
      if (std::fabs(A[r * n + c]) > best) { best = std::fabs(A[r * n + c]); piv = r; }
    }
    if (best < 1e-300) return false;
    if (piv != c) {
      for (int k = 0; k < n; ++k) { std::swap(A[c * n + k], A[piv * n + k]); std::swap(I[c * n + k], I[piv * n + k]); }
    }
    double d = 1.0 / A[c * n + c];
    for (int k = 0; k < n; ++k) { A[c * n + k] *= d; I[c * n + k] *= d; }
    for (int r = 0; r < n; ++r) {
      if (r == c) continue;
      double f = A[r * n + c];
      if (f == 0) continue;
      for (int k = 0; k < n; ++k) { A[r * n + k] -= f * A[c * n + k]; I[r * n + k] -= f * I[c * n + k]; }
    }
  }
  A.swap(I);
  return true;
}

TrackFitter::TrackFitter(const SoLIDGeometry_service* geo, const SoLIDField_service* field, double step_cm)
    : m_geo(geo), m_field(field), m_prop(field, step_cm) {}

bool TrackFitter::predict(const TrackParams& par, const std::vector<FitHit>& hits, std::vector<TrackState>& states) const {
  states.clear();
  TrackState s;
  s.pos = {0, 0, par.zv};
  s.mom = par.momentum();
  s.charge = par.charge();
  for (const auto& h : hits) {
    if (!m_prop.propagateToZ(s, h.z)) return false;
    states.push_back(s);
  }
  return true;
}

bool TrackFitter::seedFromTwoHits(const FitHit& a, const FitHit& b, double b_eff, TrackParams& out) const {
  // vertex z from the (R,z) straight line
  const double ra = std::hypot(a.x, a.y), rb = std::hypot(b.x, b.y);
  if (std::fabs(rb - ra) < 1e-6) return false;
  double zv = a.z - ra * (b.z - a.z) / (rb - ra);
  // circle through origin, a and b: 2 xa xc + 2 ya yc = ra^2 ; 2 xb xc + 2 yb yc = rb^2
  const double det = 4 * (a.x * b.y - a.y * b.x);
  if (std::fabs(det) < 1e-9) {
    // collinear with the origin -> straight track (very high momentum)
    out.qoverp = -1e-3;
    out.theta = std::atan2(rb, b.z - zv);
    out.phi = std::atan2(b.y, b.x);
    out.zv = zv;
    return true;
  }
  const double xc = (ra * ra * 2 * b.y - rb * rb * 2 * a.y) / det;
  const double yc = (2 * a.x * rb * rb - 2 * b.x * ra * ra) / det;
  const double Rc = std::hypot(xc, yc);
  const double pT = 0.00299792458 * b_eff * Rc;  // GeV
  // sense of rotation: origin -> a -> b
  const double cross = a.x * (b.y - a.y) - a.y * (b.x - a.x);
  const int charge = (cross > 0) ? -1 : +1;  // CCW in Bz>0 <=> negative charge
  // arc length from the origin to b
  const double chord = rb;
  double half = chord / (2 * Rc);
  if (half > 1) half = 1;
  const double arc = 2 * Rc * std::asin(half);
  // polar angle from the straight line vertex -> first point (the solenoid bends in phi); the arc
  // length is only used as a cross-check for very low momenta
  double theta = std::atan2(ra, a.z - zv);
  if (!(theta > 0)) theta = std::atan2(arc, b.z - zv);
  // tangent at the origin
  Vec3 t1(yc, -xc, 0), t2(-yc, xc, 0);
  Vec3 ta(a.x, a.y, 0);
  Vec3 t = (t1.dot(ta) > t2.dot(ta)) ? t1 : t2;
  out.phi = std::atan2(t.y, t.x);
  out.theta = theta;
  out.zv = zv;
  const double p = pT / std::sin(theta);
  if (!(p > 0.05) || !std::isfinite(p)) return false;
  out.qoverp = charge / p;
  return true;
}

void TrackFitter::buildScatterers(const TrackParams& par, const std::vector<FitHit>& hits,
                                  const std::vector<TrackState>& states, std::vector<Scatterer>& out) const {
  out.clear();
  const double p = par.p();
  const auto& tg = m_geo->target();
  Vec3 u0 = par.momentum().unit();
  const double sinth = std::max(1e-6, std::sin(par.theta)), costh = std::max(1e-6, std::cos(par.theta));

  // --- target: LH2 path until the cell wall or the exit window, plus the Al wall/window
  double L_side = tg.radius / sinth;
  double L_end = (tg.z_center + tg.half_length - par.zv) / costh;
  double x_over_x0;
  double L_t;
  if (L_side < L_end) {
    L_t = L_side;
    x_over_x0 = L_side / tg.x0_lh2 + tg.wall_thickness / (tg.x0_al * sinth);
  } else {
    L_t = L_end;
    x_over_x0 = L_end / tg.x0_lh2 + tg.window_thickness / (tg.x0_al * costh);
  }
  double z_exit = par.zv + L_t * costh;
  out.push_back({z_exit, L_t, highlandTheta0(p, 1.0, x_over_x0), u0});

  // --- air between the target exit and the first plane, and between planes: two Gauss points per gap
  double z_prev = z_exit, path_prev = L_t;
  Vec3 u_prev = u0;
  const double air_x0 = m_geo->airX0();
  for (size_t i = 0; i < hits.size(); ++i) {
    const auto& st = states[i];
    const double dz = st.pos.z - z_prev;
    const double dpath = st.path - path_prev;
    if (dz > 0.5) {
      const double f1 = 0.5 - 0.5 / std::sqrt(3.0), f2 = 0.5 + 0.5 / std::sqrt(3.0);
      const double th = highlandTheta0(p, 1.0, 0.5 * dpath / air_x0);
      out.push_back({z_prev + f1 * dz, path_prev + f1 * dpath, th, u_prev});
      out.push_back({z_prev + f2 * dz, path_prev + f2 * dpath, th, u_prev});
    }
    // the GEM plane itself (only affects downstream planes); calorimeter points carry no material
    Vec3 u = st.dir();
    if (hits[i].plane > 0) {
      const auto& g = m_geo->gemPlane(hits[i].plane);
      out.push_back({st.pos.z + 1e-3, st.path, highlandTheta0(p, 1.0, g.x0_fraction / std::max(1e-3, u.z)), u});
    }
    z_prev = st.pos.z;
    path_prev = st.path;
    u_prev = u;
  }
}

void TrackFitter::measurementCovariance(const TrackParams& par, const std::vector<FitHit>& hits,
                                        const std::vector<TrackState>& states, std::vector<double>& V) const {
  const int N = static_cast<int>(hits.size());
  const int n = 2 * N;
  V.assign(n * n, 0.0);
  for (int i = 0; i < N; ++i) {
    V[(2 * i) * n + 2 * i] = hits[i].sx * hits[i].sx;
    V[(2 * i + 1) * n + 2 * i + 1] = hits[i].sy * hits[i].sy;
  }
  std::vector<Scatterer> sc;
  buildScatterers(par, hits, states, sc);
  // kink at scatterer s -> displacement at plane j: (dx,dy) = M_sj (a1,a2)
  std::vector<double> Mx1(N), Mx2(N), My1(N), My2(N);
  for (const auto& s : sc) {
    if (s.theta0 <= 0) continue;
    Vec3 e1, e2;
    TrackPropagator::transverseBasis(s.u, e1, e2);
    const double uz = std::max(1e-3, s.u.z);
    Vec3 d1 = e1 - s.u * (e1.z / uz);  // shift direction in the plane for a kink along e1
    Vec3 d2 = e2 - s.u * (e2.z / uz);
    for (int j = 0; j < N; ++j) {
      if (states[j].pos.z <= s.z) { Mx1[j] = Mx2[j] = My1[j] = My2[j] = 0; continue; }
      const double L = states[j].path - s.path;
      Mx1[j] = L * d1.x; My1[j] = L * d1.y;
      Mx2[j] = L * d2.x; My2[j] = L * d2.y;
    }
    const double t2 = s.theta0 * s.theta0;
    for (int j = 0; j < N; ++j) {
      for (int k = 0; k < N; ++k) {
        V[(2 * j) * n + 2 * k] += t2 * (Mx1[j] * Mx1[k] + Mx2[j] * Mx2[k]);
        V[(2 * j) * n + 2 * k + 1] += t2 * (Mx1[j] * My1[k] + Mx2[j] * My2[k]);
        V[(2 * j + 1) * n + 2 * k] += t2 * (My1[j] * Mx1[k] + My2[j] * Mx2[k]);
        V[(2 * j + 1) * n + 2 * k + 1] += t2 * (My1[j] * My1[k] + My2[j] * My2[k]);
      }
    }
  }
}

FitResult TrackFitter::fit(const TrackParams& seed, const std::vector<FitHit>& hits, double zv_prior,
                           double zv_prior_sigma, int max_iter) const {
  FitResult res;
  const int N = static_cast<int>(hits.size());
  if (N < 2) return res;
  const bool use_prior = zv_prior_sigma > 0;
  const int n = 2 * N + (use_prior ? 1 : 0);
  TrackParams par = seed;
  std::vector<TrackState> states;
  if (!predict(par, hits, states)) return res;

  std::vector<double> W;  // weight matrix = V^-1
  auto computeWeights = [&]() {
    std::vector<double> V;
    measurementCovariance(par, hits, states, V);
    if (use_prior) {
      // extend by one row/column for the prior
      std::vector<double> Vx(n * n, 0.0);
      for (int i = 0; i < 2 * N; ++i)
        for (int j = 0; j < 2 * N; ++j) Vx[i * n + j] = V[i * (2 * N) + j];
      Vx[(n - 1) * n + n - 1] = zv_prior_sigma * zv_prior_sigma;
      V.swap(Vx);
    }
    if (!invertMatrix(V, n)) return false;
    W.swap(V);
    return true;
  };
  if (!computeWeights()) return res;

  auto residuals = [&](const TrackParams& p, const std::vector<TrackState>& st, std::vector<double>& r) {
    r.resize(n);
    for (int i = 0; i < N; ++i) {
      r[2 * i] = hits[i].x - st[i].pos.x;
      r[2 * i + 1] = hits[i].y - st[i].pos.y;
    }
    if (use_prior) r[n - 1] = zv_prior - p.zv;
  };
  auto chi2of = [&](const std::vector<double>& r) {
    double c = 0;
    for (int i = 0; i < n; ++i) {
      double s = 0;
      for (int j = 0; j < n; ++j) s += W[i * n + j] * r[j];
      c += r[i] * s;
    }
    return c;
  };

  std::vector<double> r;
  residuals(par, states, r);
  double chi2 = chi2of(r);
  const double dpar[4] = {1e-4, 1e-4, 1e-4, 0.1};  // finite-difference steps for q/p, theta, phi, zv
  std::vector<double> J(n * 4);
  int iter = 0;
  for (iter = 0; iter < max_iter; ++iter) {
    // Jacobian of the residuals w.r.t. the parameters
    for (int k = 0; k < 4; ++k) {
      TrackParams pp = par;
      double h = dpar[k];
      if (k == 0) h = std::max(1e-5, 1e-3 * std::fabs(par.qoverp));
      (k == 0 ? pp.qoverp : k == 1 ? pp.theta : k == 2 ? pp.phi : pp.zv) += h;
      std::vector<TrackState> stp;
      if (!predict(pp, hits, stp)) return res;
      std::vector<double> rp;
      residuals(pp, stp, rp);
      for (int i = 0; i < n; ++i) J[i * 4 + k] = (rp[i] - r[i]) / h;  // d r / d par
    }
    // normal equations: (J^T W J) delta = -J^T W r
    std::vector<double> A(16, 0.0), b(4, 0.0);
    std::vector<double> WJ(n * 4, 0.0), Wr(n, 0.0);
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        const double w = W[i * n + j];
        if (w == 0) continue;
        Wr[i] += w * r[j];
        for (int k = 0; k < 4; ++k) WJ[i * 4 + k] += w * J[j * 4 + k];
      }
    }
    for (int k = 0; k < 4; ++k) {
      for (int l = 0; l < 4; ++l) {
        double s = 0;
        for (int i = 0; i < n; ++i) s += J[i * 4 + k] * WJ[i * 4 + l];
        A[k * 4 + l] = s;
      }
      double s = 0;
      for (int i = 0; i < n; ++i) s += J[i * 4 + k] * Wr[i];
      b[k] = -s;
    }
    std::vector<double> Ainv = A;
    if (!invertMatrix(Ainv, 4)) return res;
    double delta[4];
    for (int k = 0; k < 4; ++k) {
      delta[k] = 0;
      for (int l = 0; l < 4; ++l) delta[k] += Ainv[k * 4 + l] * b[l];
    }
    // damped update with step halving if chi2 increases
    double lambda = 1.0;
    bool improved = false;
    for (int attempt = 0; attempt < 6; ++attempt) {
      TrackParams pn = par;
      pn.qoverp += lambda * delta[0];
      pn.theta += lambda * delta[1];
      pn.phi += lambda * delta[2];
      pn.zv += lambda * delta[3];
      if (pn.theta <= 0.01 || pn.theta >= 1.5 || std::fabs(pn.qoverp) > 20) { lambda *= 0.5; continue; }
      std::vector<TrackState> stn;
      if (!predict(pn, hits, stn)) { lambda *= 0.5; continue; }
      std::vector<double> rn;
      residuals(pn, stn, rn);
      double c = chi2of(rn);
      if (c <= chi2 + 1e-9) {
        double dchi2 = chi2 - c;
        par = pn; states.swap(stn); r.swap(rn); chi2 = c;
        improved = true;
        if (iter == 1 || iter == 3) { if (!computeWeights()) return res; chi2 = chi2of(r); }
        if (dchi2 < 1e-3 && iter > 1) { iter = max_iter; }
        break;
      }
      lambda *= 0.5;
    }
    if (!improved) break;
    for (int k = 0; k < 16; ++k) res.cov[k] = Ainv[k];
  }
  if (res.cov[0] == 0) {
    // compute covariance at the final point if the loop ended before storing it
    std::vector<double> A(16, 0.0);
    for (int k = 0; k < 4; ++k)
      for (int l = 0; l < 4; ++l) {
        double s = 0;
        for (int i = 0; i < n; ++i)
          for (int j = 0; j < n; ++j) s += J[i * 4 + k] * W[i * n + j] * J[j * 4 + l];
        A[k * 4 + l] = s;
      }
    if (!invertMatrix(A, 4)) return res;
    for (int k = 0; k < 16; ++k) res.cov[k] = A[k];
  }
  res.ok = std::isfinite(chi2);
  res.params = par;
  res.chi2 = chi2;
  res.ndf = n - 4;
  res.iterations = iter;
  res.pulls.resize(n);
  for (int i = 0; i < n; ++i) {
    double var = (W[i * n + i] > 0) ? 1.0 / W[i * n + i] : 1.0;
    res.pulls[i] = r[i] / std::sqrt(var);
  }
  return res;
}

}  // namespace solid
