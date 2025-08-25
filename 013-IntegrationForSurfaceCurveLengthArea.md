# 📘 1. 배경 – 적분과 ODE의 차이

### 적분 문제

$\[ I = \int_a^b f(x)\, dx \]$

- 함수 \(f(x)\)를 여러 점에서 평가해 면적을 근사.


### 초기값 문제 (ODE)

$\[ y'(t) = f(t,y), \quad y(t_0) = y_0 \]$

- 적분과 비슷해 보이지만, 여기서는 미지수 $\(y(t)\)$ 자체를 수치적으로 추적해야 함.  
- 즉, $\(f\)$ 를 직접 적분하는 게 아니라, **미분방정식의 해를 적분적 접근으로 근사**.

---

# 📘 2. Trapezoidal / Simpson과 ODE 해법

### Trapezoidal rule applied to ODE

$\[ y_{n+1} = y_n + \frac{h}{2} \big( f(t_n,y_n) + f(t_{n+1}, y_{n+1}) \big) \]$

- Implicit method: \(y_{n+1}\)이 양변에 있음 → 해 찾기 위해 보통 Newton iteration 필요.

---

### Simpson’s rule applied to ODE

$\[ y_{n+1} = y_n + \frac{h}{6} \Big(f(t_n,y_n) + 4 f(t_{n+\tfrac{1}{2}}, y_{n+\tfrac{1}{2}}) + f(t_{n+1}, y_{n+1}) \Big) \]$

- 역시 **암시적(implicit)** → $\(y_{n+1}\)$ 을 포함하므로 계산이 번거롭다.

---

# 📘 3. Runge–Kutta (Explicit RK)

Runge–Kutta는 사실상 Simpson/Trapezoid 아이디어를 ODE에 맞게 **명시적(explicit)**으로 재구성한 것.

### 4차 Runge–Kutta (RK4) 수식

스텝 크기 $\(h\)$ , 현재 상태 $\(y_n\)$ , 시간 $\(t_n\)$ :

$$
\[
\begin{aligned}
k_1 &= f(t_n, y_n), \\
k_2 &= f\left(t_n + \tfrac{h}{2}, \, y_n + \tfrac{h}{2}k_1\right), \\
k_3 &= f\left(t_n + \tfrac{h}{2}, \, y_n + \tfrac{h}{2}k_2\right), \\
k_4 &= f\left(t_n + h, \, y_n + h k_3\right), \\
y_{n+1} &= y_n + \tfrac{h}{6}(k_1 + 2k_2 + 2k_3 + k_4).
\end{aligned}
\]
$$

- $\(k_1, k_2, k_3, k_4\)$ 는 서로 다른 시점에서의 “슬로프(derivative)”  
- 가중합을 통해 Simpson 공식과 비슷한 포물선 보간을 흉내냄.  
- 하지만 **명시적(implicit equation 없이 계산 가능)**.

---

# 📘 4. 차이점 요약

| 구분 | Trapezoid / Simpson (적분) | Runge–Kutta (ODE) |
|------|-----------------------------|--------------------|
| 문제 | $\(\int f(x)dx\)$ | $\(y'(t) = f(t,y)\)$ |
| 근사 | 사다리꼴(1차), 포물선(2차) | 다중 슬로프 평가 후 가중합 |
| 성질 | 수치적분 전용 | ODE 해법 (적분식 재구성) |
| 암시성 | Trapezoid, Simpson은 implicit 형태로 확장됨 | RK는 explicit 방식으로 재구성 가능 |
| 정확도 | Trapezoid: 2차, Simpson: 4차 | RK4: 4차, RK45: 5차 적응적 스텝 |

---

# 📘 5. 직관 그림

- **Trapezoid**: 처음/끝 기울기 평균으로 적분  
- **Simpson**: 처음/중간/끝 세 점으로 포물선 적분  
- **RK4**: “처음, 중간, 중간, 끝” 4번의 slope 평가 → Simpson을 explicit하게 바꾼 버전

---

✅ **정리:** 

- Trapezoid/Simpson은 적분 공식을 ODE에 억지로 적용하면 implicit → 계산이 무겁다.  
- Runge–Kutta는 이를 explicit 형태로 바꿔 ODE에 맞게 설계된 방법.  
- 그래서 ODE 해석에는 RK 계열이 주로 쓰인다.

# 소스 정리
```cpp
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


  using Function1D = std::function<double(double)>;

  // ============================================================================
  // 고정 스텝 Runge-Kutta 4차 (RK4) 적분
  // ============================================================================
  static double Integrate1D_RK4(Function1D f, double a, double b, int N);

  // ============================================================================
  // Adaptive Runge-Kutta 4(5) (RK45, Dormand-Prince)
  // ============================================================================
  static double Integrate1D_RK45(Function1D f, double a, double b,
    double tol = 1e-6, int maxSteps = 10000);


  static double Integrate1D_RK4(ON_IntegrationFunction& f, void* ctx, double a, double b, int N = 1000);
  static bool   Integrate1D_RK45(ON_IntegrationFunction& f, void* ctx, double a, double b, double& out,
    double rel_tol = 1e-6, double abs_tol = 1e-9,
    double h_init = 1e-2, double h_min = 1e-12, int max_steps = 200000);
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
  tol = ON_ClampTol(tol);
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
  if (!ON_IS_FINITE(S0)) return 0.0;
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
  tol = ON_ClampTol(tol);
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
  if (!ON_IS_FINITE(S0)) return 0.0;
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

double ON_Integrator::Integrate1D_RK4(Function1D f, double a, double b, int N)
{
  double h = (b - a) / static_cast<double>(N);
  double x = a;
  double result = 0.0;

  for (int i = 0; i < N; ++i)
  {
    double k1 = f(x);
    double k2 = f(x + 0.5 * h);
    double k3 = f(x + 0.5 * h);
    double k4 = f(x + h);

    result += (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);
    x += h;
  }

  return result;
}


double ON_Integrator::Integrate1D_RK45(Function1D f, double a, double b,
  double tol, int maxSteps)
{
  double result = 0.0;
  double x = a;
  double h = (b - a) * 0.1; // 초기 스텝

  while (x < b && maxSteps-- > 0)
  {
    if (x + h > b)
      h = b - x;

    // Dormand-Prince 계수
    double k1 = f(x);
    double k2 = f(x + 0.25 * h);
    double k3 = f(x + 3.0 / 8.0 * h);
    double k4 = f(x + 12.0 / 13.0 * h);
    double k5 = f(x + h);
    double k6 = f(x + 0.5 * h);

    // 4차 및 5차 근사 (단순화된 샘플 구현)
    double y4 = h * (25.0 / 216.0 * k1 + 1408.0 / 2565.0 * k3
      + 2197.0 / 4104.0 * k4 - 1.0 / 5.0 * k5);
    double y5 = h * (16.0 / 135.0 * k1 + 6656.0 / 12825.0 * k3
      + 28561.0 / 56430.0 * k4 - 9.0 / 50.0 * k5
      + 2.0 / 55.0 * k6);

    double err = std::fabs(y5 - y4);

    if (err <= tol) // accept
    {
      result += y5;
      x += h;
    }

    // step size 조정
    double safety = 0.9;
    double power = 0.2; // 1/(order+1) = 1/5
    if (err > 0)
      h *= safety * std::pow(tol / err, power);
    else
      h *= 2.0;
  }

  return result;
}


double ON_Integrator::Integrate1D_RK4(
  ON_IntegrationFunction& f,
  void* ctx,
  double a, double b,
  int N /*=1000*/)
{
  if (a == b) return 0.0;
  if (N < 1)  N = 1;

  // 방향(부호) 처리: b < a면 뒤집어서 적분 후 부호 보정
  double sign = 1.0;
  if (b < a) { std::swap(a, b); sign = -1.0; }

  const double h = (b - a) / static_cast<double>(N);
  double x = a;
  double y = 0.0; // y(a)=0 ⇒ y(b)=∫ f

  for (int i = 0; i < N; ++i)
  {
    const double k1 = f(x, ctx);
    const double k2 = f(x + 0.5 * h, ctx);
    const double k3 = f(x + 0.5 * h, ctx);
    const double k4 = f(x + h, ctx);

    y += (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
    x += h;
  }

  return sign * y;
}

bool ON_Integrator::Integrate1D_RK45(
  ON_IntegrationFunction& f,
  void* ctx,
  double a, double b,
  double& out,
  double rel_tol /*=1e-6*/,
  double abs_tol /*=1e-9*/,
  double h_init  /*=1e-2*/,
  double h_min   /*=1e-12*/,
  int    max_steps /*=200000*/)
{
  out = 0.0;
  if (a == b) return true;

  // 방향(부호) 처리
  double sign = 1.0;
  if (b < a) { std::swap(a, b); sign = -1.0; }

  // 초기 스텝
  double h = h_init;
  if (h <= 0.0) h = (b - a) * 0.1;
  if (h > (b - a)) h = (b - a);

  // 안전계수 및 제한
  const double safety = 0.9;
  const double min_scale = 0.2;
  const double max_scale = 5.0;
  const double pow_ = 1.0 / 5.0; // 5th order controller

  double x = a;
  double y = 0.0; // y(a)=0

  int steps = 0;
  while (x < b && steps++ < max_steps)
  {
    // 마지막 스텝 조정
    if (x + h > b) h = b - x;
    if (h < h_min) { h = h_min; }

    // Dormand–Prince 5(4) tableau (y' = f(x)):
    const double k1 = f(x, ctx);
    const double k2 = f(x + (1.0 / 5.0) * h, ctx);
    const double k3 = f(x + (3.0 / 10.0) * h, ctx);
    const double k4 = f(x + (4.0 / 5.0) * h, ctx);
    const double k5 = f(x + (8.0 / 9.0) * h, ctx);
    const double k6 = f(x + h, ctx);

    // 5th-order increment (b5)
    const double incr5 =
      h * ((35.0 / 384.0) * k1
        + (500.0 / 1113.0) * k3
        + (125.0 / 192.0) * k4
        - (2187.0 / 6784.0) * k5
        + (11.0 / 84.0) * k6);

    // 4th-order increment (b4)
    const double incr4 =
      h * ((5179.0 / 57600.0) * k1
        + (7571.0 / 16695.0) * k3
        + (393.0 / 640.0) * k4
        - (92097.0 / 339200.0) * k5
        + (187.0 / 2100.0) * k6
        + (1.0 / 40.0) * f(x + h, ctx)); // same as k6 위치에서의 보정항

    const double y5 = y + incr5;
    const double y4 = y + incr4;

    // 에러 추정 및 허용오차
    const double err = std::fabs(y5 - y4);
    const double tol = (std::max)(abs_tol, rel_tol * (std::max)(std::fabs(y), 1.0));

    if (err <= tol) {
      // accept step
      y = y5;
      x += h;

      // 다음 스텝 증/감
      double factor = (err > 0.0)
        ? safety * std::pow(tol / err, pow_)
        : max_scale;
      factor = (std::min)((std::max)(factor, min_scale), max_scale);
      h *= factor;
    }
    else {
      // reject step → 스텝 줄이기
      double factor = safety * std::pow(tol / (std::max)(err, 1e-300), pow_);
      factor = (std::min)((std::max)(factor, min_scale), 1.0);
      double h_new = h * factor;
      if (h_new < h_min) {
        // 더 줄일 수 없음 → 실패
        return false;
      }
      h = h_new;
      continue; // re-try with smaller step
    }
  }

  if (x < b) {
    // max_steps 초과
    return false;
  }

  out = sign * y;
  return true;
}
```

# 샘플 정리 1
```cpp
ON_Integrator::Function1D f1 = [](double x) { return std::sin(x); };

double rk4_sin = ON_Integrator::Integrate1D_RK4(f1, 0.0, ON_PI, 1000);
double rk45_sin = ON_Integrator::Integrate1D_RK45(f1, 0.0, ON_PI, 1e-8);

std::cout << "∫ sin(x) dx [0,π] = 2" << std::endl;
std::cout << "  RK4  result = " << rk4_sin << std::endl;
std::cout << "  RK45 result = " << rk45_sin << std::endl;

// 테스트 함수 2: f(x) = exp(-x^2), 적분 [0,1] ~ 0.746824
ON_Integrator::Function1D f2 = [](double x) { return std::exp(-x * x); };

double rk4_exp = ON_Integrator::Integrate1D_RK4(f2, 0.0, 1.0, 1000);
double rk45_exp = ON_Integrator::Integrate1D_RK45(f2, 0.0, 1.0, 1e-8);

std::cout << "\n∫ exp(-x^2) dx [0,1] ≈ 0.746824" << std::endl;
std::cout << "  RK4  result = " << rk4_exp << std::endl;
std::cout << "  RK45 result = " << rk45_exp << std::endl;

```

# 결과 1
```
∫ sin(x) dx [0,π] = 2
  RK4  result = 2
  RK45 result = 2

∫ exp(-x^2) dx [0,1] ≈ 0.746824
  RK4  result = 0.746824
  RK45 result = 0.746824
```

# 샘플 정리 2
```cpp
// ----------------- 케이스 1: sin(x) [0, π] -----------------
{
  SinFunc f;
  const double a = 0.0;
  const double b = ON_PI;   // opennurbs의 PI 사용
  const double exact = 2.0;

  // RK4 (고정 스텝)
  double rk4 = ON_Integrator::Integrate1D_RK4(f, nullptr, a, b, 2000); // 스텝 수 조절 가능

  // RK45 (적응형)
  double rk45 = 0.0;
  bool ok = ON_Integrator::Integrate1D_RK45(f, nullptr, a, b, rk45, 1e-9, 1e-12, 1e-3);
  if (!ok) {
    std::printf("[sin] RK45 failed\n");
  }

  std::printf("∫ sin(x) dx, [0,π] = 2\n");
  std::printf("  RK4 : %.12f  | abs=%.3e rel=%.3e\n", rk4, abs_err(rk4, exact), rel_err(rk4, exact));
  std::printf("  RK45: %.12f  | abs=%.3e rel=%.3e\n", rk45, abs_err(rk45, exact), rel_err(rk45, exact));
}

// ----------------- 케이스 2: exp(-x^2) [0, 1] -----------------
{
  Gaussian01 f;
  const double a = 0.0;
  const double b = 1.0;
  const double reference = 0.746824132812; // 알려진 값

  double rk4 = ON_Integrator::Integrate1D_RK4(f, nullptr, a, b, 2000);
  double rk45 = 0.0;
  bool ok = ON_Integrator::Integrate1D_RK45(f, nullptr, a, b, rk45, 1e-9, 1e-12, 1e-3);
  if (!ok) {
    std::printf("[exp(-x^2)] RK45 failed\n");
  }

  std::printf("\n∫ exp(-x^2) dx, [0,1] ≈ %.12f\n", reference);
  std::printf("  RK4 : %.12f  | abs=%.3e rel=%.3e\n", rk4, abs_err(rk4, reference), rel_err(rk4, reference));
  std::printf("  RK45: %.12f  | abs=%.3e rel=%.3e\n", rk45, abs_err(rk45, reference), rel_err(rk45, reference));
}

// ----------------- 케이스 3: 파라미터 가우시안 -----------------
{
  GaussianParam f;
  struct { double c, sigma; } ctx{ 0.5, 0.10 }; // 중심 0.5, 표준편차 0.1
  const double a = 0.0;
  const double b = 1.0;

  // 이 경우 정확값을 간단히 쓰기 어려우니 RK45를 더 엄격하게 돌리면서 비교
  double rk4 = ON_Integrator::Integrate1D_RK4(f, &ctx, a, b, 4000);
  double rk45 = 0.0;
  bool ok = ON_Integrator::Integrate1D_RK45(f, &ctx, a, b, rk45, 1e-10, 1e-12, 1e-3);
  if (!ok) {
    std::printf("[GaussianParam] RK45 failed\n");
  }

  std::printf("\n∫ exp(-(x-%.2f)^2/(2*%.3f^2)) dx, [0,1]\n", ctx.c, ctx.sigma);
  std::printf("  RK4 : %.12f\n", rk4);
  std::printf("  RK45: %.12f\n", rk45);
  std::printf("  |RK45 - RK4| = %.3e\n", std::fabs(rk45 - rk4));
}


```

# 결과 정리 2

```
∫ sin(x) dx, [0,π] = 2
  RK4 : 2.000000000000  | abs=9.592e-14 rel=4.796e-14
  RK45: 2.000000000005  | abs=4.849e-12 rel=2.425e-12

∫ exp(-x^2) dx, [0,1] ≈ 0.746824132812
  RK4 : 0.746824132812  | abs=4.380e-13 rel=5.865e-13
  RK45: 0.746824133160  | abs=3.481e-10 rel=4.662e-10

∫ exp(-(x-0.50)^2/(2*0.100^2)) dx, [0,1]
  RK4 : 0.250662683757
  RK45: 0.250662683760
  |RK45 - RK4| = 2.801e-12
```


# 샘플 정리 3

```cpp

static double CalculateLengthBySimpson(FirstDerivativeLengthFunction function, 
  const ON_NurbsCurve& curve, 
  double start, 
  double end, 
  double simpson, 
  double tolearance)
{
  double length = 0.0;
  double m = (start + end) / 2.0;
  double left = ON_Integrator::Simpson(function, (void*)&curve, start, m);
  double right = ON_Integrator::Simpson(function, (void*)&curve, m, end);

  double differ = left + right - simpson;
  if (ON_AreEqual(differ, 0.0) || ON_AreLessThan(abs(differ) / 10.0, tolearance))
  {
    length = left + right + differ / 10.0;
  }
  else
  {
    length = CalculateLengthBySimpson(function, curve, start, m, left, tolearance / 2.0) 
      + CalculateLengthBySimpson(function, curve, m, end, right, tolearance / 2.0);
  }
  return length;
}

double ON_CurveApproximateLength(const ON_NurbsCurve& curve, IntegratorType type)
{
  if (curve.IsLine())
  {
    ON_3dPoint startPoint = ON_3dPoint(curve.ControlPoint(0));
    ON_3dPoint endPoint = ON_3dPoint(curve.ControlPoint(curve.CVCount() - 1));
    return startPoint.DistanceTo(endPoint);
  }

  ON_NurbsCurve reCurve = curve;
  reCurve.SetDomain(0.0, 1.0);

  int degree = reCurve.Degree();
  const double* knotVector = reCurve.Knot();
  int cntKnot = reCurve.KnotCount();
 
  double length = 0.0;
  switch (type)
  {
  case IntegratorType::Simpson:
  {
    double start = knotVector[0];
    double end = knotVector[cntKnot - 1];
    FirstDerivativeLengthFunction function;
    double simpson = 0;
    std::vector<ON_NurbsCurve> curves;
    reCurve.Decompose(curves);
    for (int i = 0; i < curves.size(); i++)
    {
      simpson += ON_Integrator::Simpson(function, (void*)&curves[i], start, end);
    }
    length = CalculateLengthBySimpson(function, reCurve, start, end, simpson, ON_SQRT_EPSILON);
    break;
  }
  case IntegratorType::GaussLegendre:
  {
    std::vector<ON_NurbsCurve> curves;
    reCurve.Decompose(curves);
    for (int i = 0; i < curves.size(); i++)
    {
      const ON_NurbsCurve& subCurve = curves[i];
      const double*  aKnots = subCurve.Knot();
      int cntLocalKnot = subCurve.KnotCount();
      double a = aKnots[0];
      double b = aKnots[cntLocalKnot - 1];
      double coefficient = (b - a) / 2.0;

      double bLength = 0.0;
      auto abscissae = GaussLegendreAbscissae;
      int size = (int)abscissae.size();
      for (int j = 0; j < size; j++)
      {
        double t = coefficient * abscissae[j] + (a + b) / 2.0;
        ON_3dPoint point;
        ON_3dVector tan;
        subCurve.TangentAt(t, point, tan);
        double derLength = tan.Length();
        if (std::isnan(derLength)) derLength = 0.0;
        bLength += GaussLegendreWeights[j] * derLength;
      }
      bLength = coefficient * bLength;
      length += bLength;
    }
    break;
  }
  case IntegratorType::Chebyshev:
  {
    std::vector<double> series = ON_Integrator::ChebyshevSeries();
    int nCntCtrlPt = reCurve.CVCount();
    for (int i = degree-1; i < nCntCtrlPt; i++)
    {
      double a = knotVector[i];
      double b = knotVector[i + 1];
      FirstDerivativeLengthFunction function;
      length += ON_Integrator::ClenshawCurtisQuadrature(function, (void*)&reCurve, a, b, series);
    }
    break;
  }
  default:
    break;
  }
  return length;
}

double ON_NurbsSurfaceApproximateArea(
  const ON_NurbsSurface& surface, 
  IntegratorType type)
{
  ON_NurbsSurface retSurface = surface;
  retSurface.SetDomain(0, 0.0, 1.0);
  retSurface.SetDomain(1, 0.0, 1.0);

  int degreeU = retSurface.Degree(0);
  int degreeV = retSurface.Degree(1);

  const double* knotVectorU = retSurface.Knot(0);
  const double* knotVectorV = retSurface.Knot(1);

  double area = 0.0;
  switch (type)
  {
  case IntegratorType::Simpson:
  {
    struct fun
      : ON_BinaryIntegrationFunction
    {
      double operator()(double u, double v, void* calFunc)const override
      {
        auto surface = (ON_NurbsSurface*)calFunc;
        ON_3dPoint S;
        ON_3dVector Su, Sv;
        surface->TangentAt(u, v, S, Su, Sv);
        double E = ON_DotProduct(Su, Su);
        double F = ON_DotProduct(Su, Sv);
        double G = ON_DotProduct(Sv, Sv);
        double ds = std::sqrt(E * G - F * F);
        return ds;
      }
    } function;

    // the initial search range
    struct UVRange
    {
      double u1, u2, v1, v2;
      double area;
    };

    UVRange init;
    init.u1 = 0;
    init.u2 = 1;
    init.v1 = 0;
    init.v2 = 1;
    init.area = ON_Integrator::Simpson(function, (void*)&retSurface, init.u1, init.u2, init.v1, init.v2);
    std::vector<UVRange> stack(1, init);

    while (!stack.empty())
    {
      UVRange uva = stack.back();
      stack.pop_back();
      // Bisect uv into 4 parts.
      double du = uva.u2 - uva.u1;
      double dv = uva.v2 - uva.v1;
      double hdu = 0.5 * du;
      double hdv = 0.5 * dv;
      UVRange uva1, uva2, uva3, uva4;
      uva1.u1 = uva.u1;
      uva1.u2 = uva.u1 + hdu;
      uva1.v1 = uva.v1;
      uva1.v2 = uva.v1 + hdv;
      uva1.area = ON_Integrator::Simpson(function, (void*)&retSurface, uva1.u1, uva1.u2, uva1.v1, uva1.v2);
      uva2.u1 = uva.u1 + hdu;
      uva2.u2 = uva.u2;
      uva2.v1 = uva.v1;
      uva2.v2 = uva.v1 + hdv;
      uva2.area = ON_Integrator::Simpson(function, (void*)&retSurface, uva2.u1, uva2.u2, uva2.v1, uva2.v2);
      uva3.u1 = uva.u1;
      uva3.u2 = uva.u1 + hdu;
      uva3.v1 = uva.v1 + hdv;
      uva3.v2 = uva.v2;
      uva3.area = ON_Integrator::Simpson(function, (void*)&retSurface, uva3.u1, uva3.u2, uva3.v1, uva3.v2);
      uva4.u1 = uva.u1 + hdu;
      uva4.u2 = uva.u2;
      uva4.v1 = uva.v1 + hdv;
      uva4.v2 = uva.v2;
      uva4.area = ON_Integrator::Simpson(function, (void*)&retSurface, uva4.u1, uva4.u2, uva4.v1, uva4.v2);

      // sum area
      double areaNew = uva1.area + uva2.area + uva3.area + uva4.area;

      // The error is tolerated.
      if (std::fabs(areaNew - uva.area) < ON_SQRT_FLOAT_EPSILON)
      {
        // Accumulate to the final area.
        area += areaNew;
      }
      else
      {
        // Continue bisections.
        stack.push_back(uva1);
        stack.push_back(uva2);
        stack.push_back(uva3);
        stack.push_back(uva4);
      }
    }

    break;
  }
  case IntegratorType::GaussLegendre:
  {
    auto& abscissae = GaussLegendreAbscissae;

    ON_SimpleArray<ON_NurbsSurface*> decomposeSurfaces;
    ON_DecomposeSurface(retSurface, decomposeSurfaces);
    for (int i = 0; i < decomposeSurfaces.Count(); i++)
    {
      ON_NurbsSurface* deSurface = decomposeSurfaces[i];
      const double* aKnotsU = deSurface->Knot(0);
      const double* aKnotsV = deSurface->Knot(1);

      int cntKnotU = deSurface->KnotCount(0);
      int cntKnotV = deSurface->KnotCount(1);

      double a = aKnotsU[0];
      double b = aKnotsU[cntKnotU - 1];
      double coefficient1 = (b - a) * 0.5;

      double c = aKnotsV[0];
      double d = aKnotsV[cntKnotV - 1];
      double coefficient2 = (d - c) * 0.5;

      double bArea = 0.0;
      auto& abscissae = GaussLegendreAbscissae;
      int size = (int)abscissae.size();
      for (int k = 0; k < size; k++)
      {
        double u = coefficient1 * abscissae[k] + (a + b) * 0.5;
        for (int j = 0; j < size; j++)
        {
          double v = coefficient2 * abscissae[j] + (c + d) * 0.5;
          ON_3dPoint S;
          ON_3dVector Su, Sv;
          deSurface->TangentAt(u, v, S, Su, Sv);
          double E = ON_DotProduct(Su, Su);
          double F = ON_DotProduct(Su, Sv);
          double G = ON_DotProduct(Sv, Sv);
          double ds = std::sqrt(E * G - F * F);
          bArea += GaussLegendreWeights[k] * GaussLegendreWeights[j] * ds;
        }
      }
      bArea = coefficient1 * coefficient2 * bArea;
      area += bArea;
      delete deSurface;
    }
    break;
  }
  case IntegratorType::Chebyshev:
  {
    struct AreaData
    {
      const ON_NurbsSurface* Surface{ nullptr };
      double parmV;
      double curKnotU;
      double nextKnotU;
      double curKnotV;
      double nextKnotV;
      const std::vector<double>& Series;

      AreaData(const ON_NurbsSurface* surface, const std::vector<double>& series) :
        Surface(surface), Series(series), parmV(0.0), curKnotU(0.0), nextKnotU(1.0),
        curKnotV(0.0), nextKnotV(0.0)
      {
      }
    };
    class AreaCoreFunction : public ON_IntegrationFunction
    {
      double operator()(double parameter, void* customFunc)
      {
        AreaData* data = (AreaData*)customFunc;
        ON_3dPoint pt;
        ON_3dVector Su;
        ON_3dVector Sv;
        data->Surface->TangentAt(parameter, data->parmV, pt, Su, Sv);
        double E = ON_DotProduct(Su, Su);
        double F = ON_DotProduct(Su, Sv);
        double G = ON_DotProduct(Sv, Sv);
        double ds = std::sqrt(E * G - F * F);
        return ds;
      }
    };
    class AreaWrapperFunction : public ON_IntegrationFunction
    {
      double operator()(double parameter, void* customData)
      {
        static std::vector<double> series;
        AreaCoreFunction areaCoreFunction;
        AreaData* data = (AreaData*)customData;
        data->parmV = parameter;
        return ON_Integrator::ClenshawCurtisQuadrature2(areaCoreFunction, data, data->curKnotU, data->nextKnotU, data->Series);
      }
    };
    
    std::vector<double> series = ON_Integrator::ChebyshevSeries();
    AreaWrapperFunction function;
    AreaData data(&retSurface, series);

    for (int i = degreeU-1; i < retSurface.CVCount(0)-1; i++)
    {
      data.curKnotU = knotVectorU[i];
      data.nextKnotU = knotVectorU[i + 1];
      for (int j = degreeV-1; j < retSurface.CVCount(1)-1; j++)
      {
        data.curKnotV = knotVectorV[j];
        data.nextKnotV = knotVectorV[j + 1];
        area += ON_Integrator::ClenshawCurtisQuadrature2(function, (void*)&data, data.curKnotV, data.nextKnotV, series);
      }
    }
    break;
  }
  default:
    break;
  }
  return area;
}




```
----
