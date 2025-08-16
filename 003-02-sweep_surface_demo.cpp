// sweep_surface_demo.cpp (standalone)
#include <algorithm>
#include <cmath>
#include <vector>
#include "opennurbs_public.h"

struct ON_SweepParams {
  int   k_sections  = 36;
  int   loft_samples= 12;
  int   degree_u    = 3;
  bool  sweep_u     = true;
  double approx_E   = 0.0;
};

static bool ON_TangentAt(const ON_Curve& c, double t, ON_3dVector& T) {
  ON_3dPoint p; ON_3dVector tan;
  if (!c.Ev1Der(t, p, tan)) return false;
  if (!tan.Unitize()) return false;
  T = tan; return true;
}

static void ON_TransportFrame(const ON_3dVector& prevT, const ON_3dVector& curT,
                              ON_3dVector& X, ON_3dVector& Y, ON_3dVector& Z) {
  ON_3dVector axis = ON_CrossProduct(prevT, curT);
  double s = axis.Length();
  double c = ON_DotProduct(prevT, curT);
  if (s < 1e-9 || c > 1.0 - 1e-12) return;
  axis.Unitize();
  double angle = std::acos(std::clamp(c, -1.0, 1.0));
  ON_Xform R; R.Rotation(angle, axis, ON_3dPoint::Origin);
  X = R * X; Y = R * Y; Z = R * Z;
}

static void ON_AdaptiveParameters(const ON_Curve& c, double E, int minN, std::vector<double>& t) {
  double a, b; c.GetDomain(&a, &b);
  t.clear();
  if (E <= 0.0) {
    int N = std::max(minN, 2);
    t.resize(N);
    for (int i = 0; i < N; ++i) t[i] = a + (b - a) * (double)i / (N - 1);
    return;
  }
  struct Seg { double a, b; };
  std::vector<Seg> st = { {a, b} };
  t.reserve(256);
  while (!st.empty()) {
    auto seg = st.back(); st.pop_back();
    double m = 0.5 * (seg.a + seg.b);
    ON_3dPoint pA = c.PointAt(seg.a);
    ON_3dPoint pB = c.PointAt(seg.b);
    ON_3dPoint pM = c.PointAt(m);
    ON_3dPoint pLin = 0.5 * (pA + pB);
    if (pM.DistanceTo(pLin) > E) {
      st.push_back({seg.a, m});
      st.push_back({m, seg.b});
    } else {
      t.push_back(seg.a);
      t.push_back(seg.b);
    }
  }
  std::sort(t.begin(), t.end());
  t.erase(std::unique(t.begin(), t.end()), t.end());
  if ((int)t.size() < minN) {
    int N = std::max(minN, 2);
    t.resize(N);
    for (int i = 0; i < N; ++i) t[i] = a + (b - a) * (double)i / (N - 1);
  }
}

static bool ON_BuildFrames(const ON_Curve& curve, const std::vector<double>& tvals,
                           const ON_3dVector& Z0, std::vector<ON_Plane>& frames) {
  if (tvals.empty()) return false;
  frames.resize(tvals.size());

  ON_3dPoint O = curve.PointAt(tvals[0]);
  ON_3dVector T; if (!ON_TangentAt(curve, tvals[0], T)) return false;

  ON_3dVector Z = Z0; if (!Z.Unitize()) Z = ON_3dVector::ZAxis;
  ON_3dVector X = T;
  ON_3dVector Y = ON_CrossProduct(Z, X);
  if (!Y.Unitize()) {
    Z = (std::fabs(T.x) < 0.9) ? ON_3dVector::XAxis : ON_3dVector::YAxis;
    Y = ON_CrossProduct(Z, X); Y.Unitize();
  }
  Z = ON_CrossProduct(X, Y); Z.Unitize();
  frames[0].CreateFromFrame(O, X, Y);

  ON_3dVector prevT = T;
  for (size_t i = 1; i < tvals.size(); ++i) {
    O = curve.PointAt(tvals[i]);
    ON_3dVector curT; if (!ON_TangentAt(curve, tvals[i], curT)) return false;

    X = frames[i-1].xaxis; Y = frames[i-1].yaxis; Z = frames[i-1].zaxis;
    ON_TransportFrame(prevT, curT, X, Y, Z);

    X = curT;
    Y = Y - ON_DotProduct(Y, X) * X; Y.Unitize();
    Z = ON_CrossProduct(X, Y); Z.Unitize();
    frames[i].CreateFromFrame(O, X, Y);
    prevT = curT;
  }
  return true;
}

static ON_NurbsCurve* ON_InstanceSection(const ON_Curve& baseC, const ON_Plane& F,
                                         const ON_3dVector* scale_or_null) {
  ON_NurbsCurve* c = baseC.NurbsCurve();
  if (!c) return nullptr;

  ON_Xform M = ON_Xform::IdentityTransformation;
  ON_Xform T = ON_Xform::TranslationTransformation(F.origin - ON_3dPoint::Origin);
  M = T * M;
  ON_Xform R; R.Rotation(ON_Plane::World_xy, F);
  M = R * M;
  if (scale_or_null) {
    ON_Xform S = ON_Xform::DiagonalTransformation(ON_3dVector(
      std::max(1e-12, scale_or_null->x),
      std::max(1e-12, scale_or_null->y),
      std::max(1e-12, scale_or_null->z)));
    M = S * M;
  }
  c->Transform(M);
  return c;
}

static bool ON_LoftCurvesSimple(const std::vector<const ON_Curve*>& curves, int samples_u,
                                bool sweep_u, int deg_sweep, ON_NurbsSurface& srf) {
  if (curves.size() < 2) return false;
  double cu0, cu1; curves[0]->GetDomain(&cu0, &cu1);
  const int W = (int)curves.size();
  const int H = std::max(4, samples_u);
  std::vector<ON_3dPoint> grid(W * H);
  for (int j = 0; j < H; ++j) {
    double t = cu0 + (cu1 - cu0) * (double)j / (H - 1);
    for (int i = 0; i < W; ++i)
      grid[j * W + i] = curves[i]->PointAt(t);
  }
  int nu = sweep_u ? W : H;
  int nv = sweep_u ? H : W;
  auto at = [&](int iu, int iv) -> const ON_3dPoint& {
    return sweep_u ? grid[iv * W + iu] : grid[iu * W + iv];
  };
  const int du = std::clamp(deg_sweep, 1, nu - 1);
  const int dv = std::min(3, nv - 1);
  srf = ON_NurbsSurface(3, false, du + 1, dv + 1, nu, nv);
  for (int iv = 0; iv < nv; ++iv)
    for (int iu = 0; iu < nu; ++iu)
      srf.SetCV(iu, iv, at(iu, iv));
  srf.MakeClampedUniformKnotVector(0, 1.0);
  srf.MakeClampedUniformKnotVector(1, 1.0);
  ON_TextLog log; return srf.IsValid(&log);
}

static bool ON_CurveDomain(const ON_Curve& c, double& t0, double& t1) {
  ON_Interval d = c.Domain();
  if (!d.IsIncreasing()) return false;
  t0 = d.Min(); t1 = d.Max(); return true;
}

static ON_NurbsCurve* MakeCircle(double radius, const ON_3dPoint& center,
                                 const ON_3dVector& normal) {
  ON_Plane pln(center, normal);
  ON_Circle C(pln, radius);
  auto* c = new ON_NurbsCurve();
  if (!C.GetNurbForm(*c)) { delete c; return nullptr; }
  return c;
}

static bool ON_SweepSurfaceBasic(const ON_Curve& trajectory, const ON_Curve& section,
                                 const ON_3dVector& Z0, const ON_Curve* scaleXYZ_or_null,
                                 const ON_SweepParams& sp, ON_NurbsSurface& out_srf) {
  std::vector<double> tvals;
  if (sp.approx_E > 0.0)
    ON_AdaptiveParameters(trajectory, sp.approx_E, std::max(4, sp.k_sections + 1), tvals);
  else {
    double a, b; if (!ON_CurveDomain(trajectory, a, b)) return false;
    int K = std::max(2, sp.k_sections + 1);
    tvals.resize(K);
    for (int i = 0; i < K; ++i) tvals[i] = a + (b - a) * (double)i / (K - 1);
  }
  std::vector<ON_Plane> frames;
  if (!ON_BuildFrames(trajectory, tvals, Z0, frames)) return false;
  std::vector<const ON_Curve*> placed;
  placed.reserve(tvals.size());
  ON_SimpleArray<ON_NurbsCurve*> owned;
  for (size_t i = 0; i < tvals.size(); ++i) {
    ON_3dVector sc(1,1,1);
    if (scaleXYZ_or_null) {
      ON_3dPoint s = scaleXYZ_or_null->PointAt(tvals[i]);
      sc = ON_3dVector(s.x, s.y, s.z);
    }
    ON_NurbsCurve* c = ON_InstanceSection(section, frames[i],
                          scaleXYZ_or_null ? &sc : nullptr);
    if (!c) { for (int k = 0; k < owned.Count(); ++k) delete owned[k]; return false; }
    owned.Append(c); placed.push_back(c);
  }
  bool ok = ON_LoftCurvesSimple(placed, sp.loft_samples, sp.sweep_u, sp.degree_u, out_srf);
  for (int k = 0; k < owned.Count(); ++k) delete owned[k];
  return ok;
}

int main() {
  ON::Begin();
  ON_TextLog log;

  ON_Polyline pl;
  pl.Append(ON_3dPoint(-10,  0, 0));
  pl.Append(ON_3dPoint( -8,  3, 3));
  pl.Append(ON_3dPoint( -2,  0, 5));
  pl.Append(ON_3dPoint(  5, -3, 3));
  pl.Append(ON_3dPoint( 10,  0, 0));
  auto* traj_pl = new ON_PolylineCurve(pl);
  if (!traj_pl || !traj_pl->IsValid()) { log.Print("trajectory build failed.\n"); return 1; }

  auto* sect_nc = MakeCircle(1.2, ON_3dPoint::Origin, ON_3dVector::ZAxis);
  if (!sect_nc) { log.Print("section build failed.\n"); delete traj_pl; return 1; }

  double t0, t1; traj_pl->GetDomain(&t0, &t1);
  ON_Polyline sc_pl;
  sc_pl.Append(ON_3dPoint(1.0, 1.0, 1.0));
  sc_pl.Append(ON_3dPoint(0.4, 0.4, 1.2));
  auto* scale_crv = new ON_PolylineCurve(sc_pl);
  scale_crv->SetDomain(t0, t1);

  ON_SweepParams sp;
  sp.k_sections   = 36;
  sp.loft_samples = 12;
  sp.degree_u     = 3;
  sp.sweep_u      = true;
  sp.approx_E     = 0.0;

  ON_NurbsSurface srf;
  bool ok = ON_SweepSurfaceBasic(*traj_pl, *sect_nc, ON_3dVector::ZAxis, scale_crv, sp, srf);
  if (!ok) {
    log.Print("ON_SweepSurfaceBasic failed.\n");
    delete traj_pl; delete sect_nc; delete scale_crv;
    ON::End(); return 1;
  }

  log.Print("Sweep OK. deg(u,v)=(%d,%d), CV=(%d x %d)\n",
            srf.Degree(0), srf.Degree(1), srf.CVCount(0), srf.CVCount(1));

  //ON_Mesh mesh;
  //ON_MeshParameters mp; mp.SetGridMaxCount(5);
  //srf.CreateMesh(mp, &mesh);
  // ON_ExportBDF bdf; bdf.Run(L\"D:\\\\Temp\\\\Test1.bdf\", &mesh);

  delete traj_pl; delete sect_nc; delete scale_crv;
  ON::End();
  return 0;
}
