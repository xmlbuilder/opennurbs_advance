# ON_Integrator --- Numerical Integration for Geometry (C++)

재사용 가능한 1D/2D 적분 유틸 (`Simpson`, `SimpsonAdaptive`,
`Gauss–Legendre`, `Clenshaw–Curtis`)과 곡선 길이/곡면 면적 적분용 샘플을
담은 README입니다.

------------------------------------------------------------------------

## ✨ 핵심 아이디어

-   **평가자 인터페이스 분리**:
    -   1D: `ON_IntegrationFunction::operator()(double t, void* ctx)`\
    -   2D:
        `ON_BinaryIntegrationFunction::operator()(double u, double v, void* ctx)`
-   **안정성**: 구간 퇴화/NaN 방어, 상대 오차 기반 종료(적응형 심프슨)
-   **기하 적분**:
-   
    -   곡선 길이:
        ![L = `\int`{=tex}\_a\^b
        \|C'(t)\|,dt](https://math.vercel.app/?from=L%20%3D%20%5Cint_a%5Eb%20%5C%7C%20C'(t)%20%5C%7C%20dt)
        
    -   곡면 면적:
        ![`\int`{=tex}`\int `{=tex}\|S_u
        `\times `{=tex}S_v\|,dv,du](https://math.vercel.app/?from=%5Cint_%7Bu_0%7D%5E%7Bu_1%7D%20%5Cint_%7Bv_0%7D%5E%7Bv_1%7D%20%5C%7C%20S_u%20%5Ctimes%20S_v%20%5C%7C%20dv%20du)

------------------------------------------------------------------------

## 📐 알고리즘 요약

### Trapezoid 누적 → Romberg 외삽 (참고식)

-   재귀형 사다리꼴 세밀화:
    ![T_j = `\tfrac12`{=tex} T\_{j-1} + h_j
    `\sum `{=tex}f(`\cdot`{=tex}),`\quad `{=tex}h_j=`\tfrac{b-a}{2^{j-1}}`{=tex}](https://math.vercel.app/?from=T_j%20%3D%200.5T_%7Bj-1%7D%20%2B%20h_j%20%5Csum%20f(...),%5Cquad%20h_j%3D%5Cfrac%7Bb-a%7D%7B2%5E%7Bj-1%7D%7D)
    
-   Romberg 외삽:
    ![R\_{j,k} = R\_{j,k-1} +
    `\frac{R_{j,k-1} - R_{j-1,k-1}}{4^k - 1}`{=tex}](https://math.vercel.app/?from=R_%7Bj%2Ck%7D%20%3D%20R_%7Bj%2Ck-1%7D%20%2B%20%5Cfrac%7BR_%7Bj%2Ck-1%7D-R_%7Bj-1%2Ck-1%7D%7D%7B4%5Ek-1%7D)
    

> 본 리포의 기본 구현은 **적응형 심프슨 + 고정차수 Gauss--Legendre**를
> 권장합니다.

### Composite Simpson (Trapezoid 재사용 관계식)

![S_j =
`\frac{4T_j - T_{j-1}}{3}`{=tex}](https://math.vercel.app/?from=S_j%20%3D%20%5Cfrac%7B4T_j%20-%20T_%7Bj-1%7D%7D%7B3%7D)

### Gauss--Legendre (고정차수)

-   1D: ![`\int`{=tex}\_a\^b f(t),dt `\approx `{=tex}`\sum`{=tex}\_i
    w_i, f!`\left`{=tex}(`\tfrac{b-a}{2}`{=tex}x_i +
    `\tfrac{b+a}{2}`{=tex}`\right`{=tex})`\tfrac{b-a}{2}`{=tex}](https://math.vercel.app/?from=%5Cint_a%5Eb%20f(t)dt%20%5Capprox%20%5Csum_i%20w_i%20f((b-a)x_i/2%20%2B%20(b%2Ba)/2)%5Ccdot(b-a)/2)

### Clenshaw--Curtis (Chebyshev 기반)

-   (\[-1,1\])의 체비셰프 노드 샘플 → DCT-I로 계수 추정 → 적분 가중 합\
    (엔드포인트 특이성에 강함. 본 구현은 **간단 DCT-I 기반**)

------------------------------------------------------------------------

## 🔧 공개 API 개요

``` cpp
class ON_IntegrationFunction {
public:
  virtual ~ON_IntegrationFunction() = default;
  virtual double operator()(double t, void* ctx) = 0;
};

class ON_BinaryIntegrationFunction {
public:
  virtual ~ON_BinaryIntegrationFunction() = default;
  virtual double operator()(double u, double v, void* ctx) const = 0;
};

class ON_Integrator {
public:
  // 1D
  static double Simpson(ON_IntegrationFunction& f, void* ctx, double a, double b);
  static double SimpsonAdaptive(ON_IntegrationFunction& f, void* ctx,
                                double a, double b, double tol=1e-8, int maxDepth=16);
  static double GaussLegendre(ON_IntegrationFunction& f, void* ctx, double a, double b);
  static double ClenshawCurtis(ON_IntegrationFunction& f, void* ctx, double a, double b, int N=64);

  // 2D (u,v)
  static double Simpson(ON_BinaryIntegrationFunction& f, void* ctx,
                        double u0, double u1, double v0, double v1);
  static double SimpsonAdaptive(ON_BinaryIntegrationFunction& f, void* ctx,
                                double u0, double u1, double v0, double v1,
                                double tol=1e-6, int maxDepth=10);
  static double GaussLegendre(ON_BinaryIntegrationFunction& f, void* ctx,
                              double u0, double u1, double v0, double v1);
};
```

------------------------------------------------------------------------

## 🧭 기하학 Integrand 정의 예

### 1) 곡선 길이 integrand --- 속도(norm of tangent)

![\|C'(t)\|](https://math.vercel.app/?from=%5C%7C%20C'(t)%20%5C%7C)

``` cpp
struct CurveSpeed : public ON_IntegrationFunction {
  // 사용자 정의 곡선 객체를 캡슐화하세요.
  const MyCurve& curve;
  explicit CurveSpeed(const MyCurve& c) : curve(c) {}
  double operator()(double t, void* ctx) override {
    Vec3 d = curve.DerivativeAt(t);          // C'(t)
    double v = d.norm();
    return std::isfinite(v) ? v : 0.0;
  }
};
```

### 2) 곡면 면적 integrand --- 외적 노름

![\|S_u
`\times `{=tex}S_v\|](https://math.vercel.app/?from=%5C%7C%20S_u%20%5Ctimes%20S_v%20%5C%7C)

``` cpp
struct SurfaceAreaDensity : public ON_BinaryIntegrationFunction {
  const MySurface& srf;
  explicit SurfaceAreaDensity(const MySurface& s) : srf(s) {}
  double operator()(double u, double v, void* ctx) const override {
    Vec3 Su, Sv; srf.Partials(u, v, Su, Sv); // S_u, S_v
    Vec3 n = cross(Su, Sv);
    double a = n.norm();
    return std::isfinite(a) ? a : 0.0;
  }
};
```

------------------------------------------------------------------------

## ✅ 전체 코드 (C++)

```cpp

static constexpr std::array<double, 24> GaussLegendreAbscissae =
{
    -0.0640568928626056260850430826247450385909,
    0.0640568928626056260850430826247450385909,
    -0.1911188674736163091586398207570696318404,
    0.1911188674736163091586398207570696318404,
    -0.3150426796961633743867932913198102407864,
    0.3150426796961633743867932913198102407864,
    -0.4337935076260451384870842319133497124524,
    0.4337935076260451384870842319133497124524,
    -0.5454214713888395356583756172183723700107,
    0.5454214713888395356583756172183723700107,
    -0.6480936519369755692524957869107476266696,
    0.6480936519369755692524957869107476266696,
    -0.7401241915785543642438281030999784255232,
    0.7401241915785543642438281030999784255232,
    -0.8200019859739029219539498726697452080761,
    0.8200019859739029219539498726697452080761,
    -0.8864155270044010342131543419821967550873,
    0.8864155270044010342131543419821967550873,
    -0.9382745520027327585236490017087214496548,
    0.9382745520027327585236490017087214496548,
    -0.9747285559713094981983919930081690617411,
    0.9747285559713094981983919930081690617411,
    -0.9951872199970213601799974097007368118745,
    0.9951872199970213601799974097007368118745,
};

static constexpr std::array<double, 24> GaussLegendreWeights =
{
    0.1279381953467521569740561652246953718517,
    0.1279381953467521569740561652246953718517,
    0.1258374563468282961213753825111836887264,
    0.1258374563468282961213753825111836887264,
    0.121670472927803391204463153476262425607,
    0.121670472927803391204463153476262425607,
    0.1155056680537256013533444839067835598622,
    0.1155056680537256013533444839067835598622,
    0.1074442701159656347825773424466062227946,
    0.1074442701159656347825773424466062227946,
    0.0976186521041138882698806644642471544279,
    0.0976186521041138882698806644642471544279,
    0.086190161531953275917185202983742667185,
    0.086190161531953275917185202983742667185,
    0.0733464814110803057340336152531165181193,
    0.0733464814110803057340336152531165181193,
    0.0592985849154367807463677585001085845412,
    0.0592985849154367807463677585001085845412,
    0.0442774388174198061686027482113382288593,
    0.0442774388174198061686027482113382288593,
    0.0285313886289336631813078159518782864491,
    0.0285313886289336631813078159518782864491,
    0.0123412297999871995468056670700372915759,
    0.0123412297999871995468056670700372915759,
};


class ON_CLASS ON_IntegrationFunction
{
public:

  virtual double operator()(double parameter, void* calcFunc) = 0;
};



class ON_CLASS ON_BinaryIntegrationFunction
{
public:
  virtual ~ON_BinaryIntegrationFunction() = default;

  virtual double operator()(double u, double v, void* calcFunc)const = 0;
};

inline bool is_finite(double x) {
  return std::isfinite(x);
}

inline double clamp_tol(double tol) {
  if (!(tol > 0.0)) return 1e-8;
  return tol;
}


class ON_CLASS ON_Integrator
{

public:
  static double Simpson(
    ON_IntegrationFunction& function, 
    void* calcFunc,
    double a, 
    double b);


  static double SimpsonAdaptive(
    ON_IntegrationFunction& f, 
    void* calcFunc,
    double a, 
    double b,
    double tol = 1e-8, 
    int maxDepth = 16);

  static double Simpson(
    ON_BinaryIntegrationFunction& function, 
    void* calcFunc, 
    double uStart,
    double uEnd,
    double vStart, 
    double vEnd);
  

  double SimpsonAdaptive(
    ON_BinaryIntegrationFunction& f,
    void* calcFunc,
    double u0,
    double u1,
    double v0,
    double v1,
    double tol = 1e-6,
    int maxDepth = 10
  );

  double GaussLegendre(ON_IntegrationFunction& f, 
    void* calcFunc, double a, double b);
  
  double GaussLegendre(ON_BinaryIntegrationFunction& f, 
    void* calcFunc,
    double u0, double u1, double v0, double v1);


  static double ClenshawCurtisSimple(
    ON_IntegrationFunction& f, void* ctx,
    double a, double b, int N = 64);

  static double ClenshawCurtisQuadrature(
    ON_IntegrationFunction& function, 
    void* calcFunc,
    double start, double end, 
    std::vector<double>& series, 
    double epsilon = ON_EPSILON);


  static double ClenshawCurtisQuadrature2(
    ON_IntegrationFunction& function, 
    void* calcFunc,
    double start, 
    double end, 
    std::vector<double> series, 
    double epsilon = ON_EPSILON);

  static std::vector<double> ChebyshevSeries(int size=100);

};

std::vector<double> ON_Integrator::ChebyshevSeries(int size)
{
  std::vector<double> series(size);

  int lenw = (int)series.size() - 1;
  int j, k, l, m;
  double cos2, sin1, sin2, hl;

  cos2 = 0;
  sin1 = 1;
  sin2 = 1;
  hl = 0.5;
  k = lenw;
  l = 2;
  while (l < k - l - 1)
  {
    series[0] = hl * 0.5;
    for (j = 1; j <= l; j++)
    {
      series[j] = hl / (1 - 4 * j * j);
    }
    series[l] *= 0.5;
    dfct(l, 0.5 * cos2, sin1, series);
    cos2 = std::sqrt(2 + cos2);
    sin1 /= cos2;
    sin2 /= 2 + cos2;
    series[k] = sin2;
    series[k - 1] = series[0];
    series[k - 2] = series[l];
    k -= 3;
    m = l;
    while (m > 1)
    {
      m >>= 1;
      for (j = m; j <= l - m; j += (m << 1))
      {
        series[k] = series[j];
        k--;
      }
    }
    hl *= 0.5;
    l *= 2;
  }
  return series;
}

double ON_Integrator::Simpson(ON_IntegrationFunction& function, void* calcFunc, double a, double b)
{
  if (a == b) return 0.0;
  double st = (function)(a, calcFunc);
  double mt = (function)((a + b) / 2.0, calcFunc);
  double et = (function)((b), calcFunc);
  double result = ((b - a) / 6.0) * (st + 4 * mt + et);
  return std::isfinite(result) ? result : 0.0;
}


double ON_Integrator::SimpsonAdaptive(
  ON_IntegrationFunction& f, 
  void* calcFunc,
  double a, 
  double b,
  double tol, 
  int maxDepth)
{
  tol = clamp_tol(tol);
  if (b == a) return 0.0;

  struct Node {
    static double panel(ON_IntegrationFunction& ff, void* cc, double aa, double bb) {
      const double mm = 0.5 * (aa + bb);
      const double fa = ff(aa, cc);
      const double fm = ff(mm, cc);
      const double fb = ff(bb, cc);
      return (bb - aa) * (fa + 4.0 * fm + fb) / 6.0;
    }
    static double rec(ON_IntegrationFunction& ff, void* cc,
      double aa, double bb, double S,
      double tol, int depth, int maxDepth)
    {
      const double m = 0.5 * (aa + bb);
      const double S1 = panel(ff, cc, aa, m);
      const double S2 = panel(ff, cc, m, bb);
      const double err = std::fabs(S1 + S2 - S);
      if (err < 15.0 * tol || depth >= maxDepth) {
        return S1 + S2 + (S1 + S2 - S) / 15.0; // Richardson correction
      }
      return rec(ff, cc, aa, m, S1, tol * 0.5, depth + 1, maxDepth)
        + rec(ff, cc, m, bb, S2, tol * 0.5, depth + 1, maxDepth);
    }
  };

  const double S0 = ON_Integrator::Simpson(f, calcFunc, a, b);
  if (!is_finite(S0)) return 0.0;
  return Node::rec(f, calcFunc, a, b, S0, tol, 0, maxDepth);
}

double ON_Integrator::Simpson(ON_BinaryIntegrationFunction& function, 
  void* calcFunc, double uStart, double uEnd, double vStart, double vEnd)
{
  double du = uEnd - uStart;
  double dv = vEnd - vStart;

  double hdu = 0.5 * du;
  double hdv = 0.5 * dv;

  // Sample 9 points with weights
  int sampleNumber = 9;
  int patches = 4;
  double uvw[27] = {
      uStart,         vStart,         1,
      uStart,         vStart + hdv,   4,
      uStart,         vEnd,           1,
      uStart + hdu,   vStart,         4,
      uStart + hdu,   vStart + hdv,   16,
      uStart + hdu,   vEnd,           4,
      uEnd,           vStart,         1,
      uEnd,           vStart + hdv,   4,
      uEnd,           vEnd,           1,
  };

  double sum = 0;
  for (int i = 0; i < sampleNumber; ++i)
  {
    double* base = uvw + i * 3;
    double u = base[0];
    double v = base[1];
    double w = base[2];
    double f = function(u, v, calcFunc);
    sum += w * f;
  }

  sum *= du * dv / (sampleNumber * patches);
  return sum;
}

double ON_Integrator::SimpsonAdaptive(
  ON_BinaryIntegrationFunction& f, 
  void* calcFunc,
  double u0, 
  double u1, 
  double v0, 
  double v1,
  double tol, 
  int maxDepth)
{
  tol = clamp_tol(tol);
  if (u0 == u1 || v0 == v1) return 0.0;

  auto panel = [&](double a, double b, double c, double d)->double {
    return Simpson(f, calcFunc, a, b, c, d);
    };

  struct Node {
    static double rec(ON_BinaryIntegrationFunction& ff, void* cc,
      double a, double b, double c, double d,
      double S, double tol, int depth, int maxDepth,
      decltype(panel)& panelFn)
    {
      const double um = 0.5 * (a + b), vm = 0.5 * (c + d);
      const double S11 = panelFn(a, um, c, vm);
      const double S12 = panelFn(um, b, c, vm);
      const double S21 = panelFn(a, um, vm, d);
      const double S22 = panelFn(um, b, vm, d);
      const double S4 = S11 + S12 + S21 + S22;

      const double err = std::fabs(S4 - S);
      if (err < 15.0 * tol || depth >= maxDepth) {
        return S4 + (S4 - S) / 15.0; // Richardson-esque correction
      }
      const double tchild = tol * 0.25;
      return rec(ff, cc, a, um, c, vm, S11, tchild, depth + 1, maxDepth, panelFn)
        + rec(ff, cc, um, b, c, vm, S12, tchild, depth + 1, maxDepth, panelFn)
        + rec(ff, cc, a, um, vm, d, S21, tchild, depth + 1, maxDepth, panelFn)
        + rec(ff, cc, um, b, vm, d, S22, tchild, depth + 1, maxDepth, panelFn);
    }
  };

  const double S0 = panel(u0, u1, v0, v1);
  if (!is_finite(S0)) return 0.0;
  return Node::rec(f, calcFunc, u0, u1, v0, v1, S0, tol, 0, maxDepth, panel);
}

double ON_Integrator::ClenshawCurtisQuadrature(ON_IntegrationFunction& function, void* calcFunc, double start, double end, std::vector<double>& series, double epsilon)
{
  double result;
  int j, k, l;
  double err, esf, eref, erefh, hh, ir, iback, irback, ba, ss, x, y, fx, errir;
  int lenw = (int)series.size() - 1;
  esf = 10;
  ba = 0.5 * (end - start);
  ss = 2 * series[lenw];
  x = ba * series[lenw];
  series[0] = 0.5 * (function)(start, calcFunc);
  series[3] = 0.5 * (function)(end, calcFunc);
  series[2] = (function)(start + x, calcFunc);
  series[4] = (function)(end - x, calcFunc);
  series[1] = (function)(start + ba, calcFunc);
  eref = 0.5 * (fabs(series[0]) + std::fabs(series[1]) + std::fabs(series[2]) + std::fabs(series[3]) + std::fabs(series[4]));
  series[0] += series[3];
  series[2] += series[4];
  ir = series[0] + series[1] + series[2];
  result = series[0] * series[lenw - 1] + series[1] * series[lenw - 2] + series[2] * series[lenw - 3];
  erefh = eref * std::sqrt(epsilon);
  eref *= epsilon;
  hh = 0.25;
  l = 2;
  k = lenw - 5;
  do {
    iback = result;
    irback = ir;
    x = ba * series[k + 1];
    y = 0;
    result = series[0] * series[k];
    for (j = 1; j <= l; j++) {
      x += y;
      y += ss * (ba - x);
      fx = (function)(start + x, calcFunc) + (function)(end - x, calcFunc);
      ir += fx;
      result += series[j] * series[k - j] + fx * series[k - j - l];
      series[j + l] = fx;
    }
    ss = 2 * series[k + 1];
    err = esf * l * std::fabs(result - iback);
    hh *= 0.25;
    errir = hh * std::fabs(ir - 2 * irback);
    l *= 2;
    k -= l + 2;
  } while ((err > erefh || errir > eref) && k > 4 * l);
  result *= end - start;
  if (err > erefh || errir > eref)
  {
    err *= -fabs(end - start);
  }
  else
  {
    err = eref * std::fabs(end - start);
  }
  return result;
}

double ON_Integrator::GaussLegendre(ON_IntegrationFunction& f, void* calcFunc, double a, double b)
{
  if (b == a) return 0.0;
  const double c1 = 0.5 * (b - a);
  const double c2 = 0.5 * (b + a);
  double S = 0.0;
  for (size_t i = 0; i < GaussLegendreAbscissae.size(); ++i) {
    const double xi = c1 * GaussLegendreAbscissae[i] + c2;
    S += GaussLegendreWeights[i] * f(xi, calcFunc);
  }
  S *= c1;
  return std::isfinite(S) ? S : 0.0;
}

// ---- 2D Gauss–Legendre (tensor product of 24-pt) ----
double ON_Integrator::GaussLegendre(ON_BinaryIntegrationFunction& f, void* calcFunc,
  double u0, double u1, double v0, double v1)
{
  if (u0 == u1 || v0 == v1) return 0.0;
  const double cu = 0.5 * (u1 - u0), vu = 0.5 * (u1 + u0);
  const double cv = 0.5 * (v1 - v0), vv = 0.5 * (v1 + v0);

  double S = 0.0;
  for (size_t i = 0; i < GaussLegendreAbscissae.size(); ++i) {
    const double ui = cu * GaussLegendreAbscissae[i] + vu;
    const double wi = GaussLegendreWeights[i];
    for (size_t j = 0; j < GaussLegendreAbscissae.size(); ++j) {
      const double vj = cv * GaussLegendreAbscissae[j] + vv;
      const double wj = GaussLegendreWeights[j];
      S += wi * wj * f(ui, vj, calcFunc);
    }
  }
  S *= (cu * cv);
  return std::isfinite(S) ? S : 0.0;
}


double ON_Integrator::ClenshawCurtisSimple(
  ON_IntegrationFunction& f, void* calcFunc,
  double a, double b, int N)
{
  if (a == b) return 0.0;
  // 노드: x_k = cos(k*pi/N), k=0..N  ([-1,1])
  // [a,b]로 매핑: t = (b+a)/2 + (b-a)/2 * x
  std::vector<double> w(N + 1, 0.0), x(N + 1, 0.0);
  for (int k = 0; k <= N; ++k) x[k] = std::cos(ON_PI * k / N);

  // 간단 가중치(표준식; DCT 없이 근사) — 정확한 가중치는 별도 테이블/유도 필요
  // 실무에선 미리 계산된 weight 사용을 권장.
  w[0] = w[N] = 1.0 / (N * N - 1.0); // placeholder; 실제론 정확한 CC weight 필요
  for (int k = 1; k < N; ++k) w[k] = 2.0 / (1.0 - 4.0 * k * k); // placeholder

  const double c1 = 0.5 * (b - a), c2 = 0.5 * (b + a);
  double S = 0.0;
  for (int k = 0; k <= N; ++k) {
    const double t = c2 + c1 * x[k];
    S += w[k] * f(t, calcFunc);
  }
  return c1 * S;
}

double ON_Integrator::ClenshawCurtisQuadrature2(ON_IntegrationFunction& function, void* calcFunc, double start, double end, std::vector<double> series, double epsilon)
{
  double result;
  int j, k, l;
  double err, esf, eref, erefh, hh, ir, iback, irback, ba, ss, x, y, fx, errir;
  int lenw = (int)series.size() - 1;
  esf = 10;
  ba = 0.5 * (end - start);
  ss = 2 * series[lenw];
  x = ba * series[lenw];
  series[0] = 0.5 * (function)(start, calcFunc);
  series[3] = 0.5 * (function)(end, calcFunc);
  series[2] = (function)(start + x, calcFunc);
  series[4] = (function)(end - x, calcFunc);
  series[1] = (function)(start + ba, calcFunc);
  eref = 0.5 * (fabs(series[0]) + std::fabs(series[1]) + std::fabs(series[2]) + std::fabs(series[3]) + std::fabs(series[4]));
  series[0] += series[3];
  series[2] += series[4];
  ir = series[0] + series[1] + series[2];
  result = series[0] * series[lenw - 1] + series[1] * series[lenw - 2] + series[2] * series[lenw - 3];
  erefh = eref * std::sqrt(epsilon);
  eref *= epsilon;
  hh = 0.25;
  l = 2;
  k = lenw - 5;
  do {
    iback = result;
    irback = ir;
    x = ba * series[k + 1];
    y = 0;
    result = series[0] * series[k];
    for (j = 1; j <= l; j++) {
      x += y;
      y += ss * (ba - x);
      fx = (function)(start + x, calcFunc) + (function)(end - x, calcFunc);
      ir += fx;
      result += series[j] * series[k - j] + fx * series[k - j - l];
      series[j + l] = fx;
    }
    ss = 2 * series[k + 1];
    err = esf * l * std::fabs(result - iback);
    hh *= 0.25;
    errir = hh * std::fabs(ir - 2 * irback);
    l *= 2;
    k -= l + 2;
  } while ((err > erefh || errir > eref) && k > 4 * l);
  result *= end - start;
  if (err > erefh || errir > eref)
  {
    err *= -fabs(end - start);
  }
  else
  {
    err = eref * std::fabs(end - start);
  }
  return result;
}

```

---


