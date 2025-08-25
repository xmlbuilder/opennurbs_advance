# 🧮 ON_IntegratorEx 알고리즘 설명 (Trapezoid → Romberg / Simpson) + 기하적 적분식

이 문서는 `ON_IntegratorEx`와 그 위에서 동작하는 **곡선 길이 적분**과 **곡면 면적 적분**의 이론식을 정리합니다.

---

## 1) Trapezoid Rule (재귀적 세밀화)

기본 사다리꼴 공식:
![trap](https://math.vercel.app/?from=T_n%20%3D%20h%5CBig%5B%5Ctfrac12%20f(a)%2B%5Csum_%7Bk%3D1%7D%5E%7Bn-1%7D%20f(a%2Bk%20h)%2B%5Ctfrac12%20f(b)%5CBig%5D%2C%5Cquad%20h%3D%5Cfrac%7Bb-a%7D%7Bn%7D)

사다리꼴 값을 점점 **세밀화**하며 갱신하는 재귀형(코드의 `Trap(e)zoid`와 일치):
![trap-refine](https://math.vercel.app/?from=T_j%20%3D%200.5%20T_%7Bj-1%7D%20%2B%20h_j%20%5Csum_%7Bk%3D1%7D%5E%7B2%5E%7Bj-2%7D%7D%20f%5Cbig(a%20%2B%20(2k-1)h_j%5Cbig)%2C%5Cquad%20h_j%20%3D%20%5Cfrac%7Bb-a%7D%7B2%5E%7Bj-1%7D%7D)

---

## 2) Romberg Integration (Richardson Extrapolation)

사다리꼴 값 $\(T_j\)$ 를 바탕으로, 격자간격 $\(h\to 0\)$ 으로의 외삽을 반복합니다.  
초기값:
![romberg-base](https://math.vercel.app/?from=R_%7Bj%2C0%7D%20%3D%20T_%7B2%5Ej%7D)

외삽 재귀식(코드의 `PolynomialInterpolation` 자리—표준식 표기):
![romberg-rec](https://math.vercel.app/?from=R_%7Bj%2Ck%7D%20%3D%20R_%7Bj%2Ck-1%7D%20%2B%20%5Cfrac%7BR_%7Bj%2Ck-1%7D%20-%20R_%7Bj-1%2Ck-1%7D%7D%7B4%5Ek%20-%201%7D)

> 구현 메모  
> • 본 코드는 `PolynomialInterpolation`(Neville 보간)으로 **$\(h=0\)$** 지점의 값을 직접 추정합니다.  
> • `QRomberg()` 내부에서 보간 입력 배열을 **1‑based**로 다루므로, 포인터 오프셋을 `+1` 해 준 점이 중요합니다.

수렴 판단(상대 오차):
![stop-rel](https://math.vercel.app/?from=%7CR_%7Bj%2Ck%7D-%20R_%7Bj%2Ck-1%7D%7C%20%3C%20%5Cmathrm%7Btol%7D%5Ccdot(%7C%20R_%7Bj%2Ck%7D%20%7C%2B1))

---

## 3) Composite Simpson (Trapezoid 재사용)

사다리꼴 두 단계 $\(T_{j-1},T_j\)$ 로 Simpson 값을 얻는 표준 관계식:
![simp-from-trap](https://math.vercel.app/?from=S_j%20%3D%20%5Cfrac%7B4T_j%20-%20T_%7Bj-1%7D%7D%7B3%7D)

코드의 `QSimpsons()`는 위 식으로 $\(S_j\)$ 를 갱신하고,  
![stop-rel-s](https://math.vercel.app/?from=%7CS_j-S_%7Bj-1%7D%7C%20%3C%20%5Cmathrm%7Btol%7D%5Ccdot(%7C%20S_%7Bj-1%7D%20%7C%2B1))
형태의 상대 오차로 종료합니다. (초기 1–2 스텝은 비교 건너뛰는 가드가 안정적)

---

## 4) 곡선 길이 적분 (Curve Length)

곡선 $\(C(t)\)$ 의 길이:
![len-int](https://math.vercel.app/?from=L%20%3D%20%5Cint_a%5Eb%20%5C%7C%20C'(t)%20%5C%7C%20%5C%2C%20dt)

구현에서 integrand는 **속도 크기** $\(\|C'(t)\|\)$ .  
`ON_CalcCurveLengthByTangent()`는 도메인 전체에 대해 **Romberg**로 적분합니다.

---

## 5) 곡면 면적 적분 (Surface Area)

매개화 곡면 \(S(u,v)\)에 대해 제1기본형:
![EGF](https://math.vercel.app/?from=E%20%3D%20S_u%5Ccdot%20S_u%2C%5Cquad%20F%20%3D%20S_u%5Ccdot%20S_v%2C%5Cquad%20G%20%3D%20S_v%5Ccdot%20S_v)

면적 소요소(두 가지 동치 표현):
![area-el](https://math.vercel.app/?from=%5C%7C%20S_u%20%5Ctimes%20S_v%20%5C%7C%20%3D%20%5Csqrt%7BEG%20-%20F%5E2%7D)

총 면적(이중 적분):
![area-int](https://math.vercel.app/?from=A%20%3D%20%5Cint_%7Bu_0%7D%5E%7Bu_1%7D%20%5Cint_%7Bv_0%7D%5E%7Bv_1%7D%20%5C%7C%20S_u%20%5Ctimes%20S_v%20%5C%7C%20%5C%2C%20dv%5C%2Cdu)

구현에서는 **외적 노름** $\(\|S_u\times S_v\|\)$ 을 integrand로 사용하고,  
바깥쪽 적분( $\(u\)$ ) 루프와 안쪽 적분( $\(v\)$ iso‑curve 길이 ) 모두 **Romberg**로 수행합니다.

---

## 6) 정확도/수렴과 방어 로직

- **사다리꼴 누적**: 이전 합을 재활용해 연산량 절감  
- **Romberg/Simpson 수렴**: 상대 오차 기준 (값이 0 근처일 땐 `+1` 가드)  
- **내부/외부 오차 배분**:  
  - 면적에서 내부 iso‑curve 길이 적분 허용오차는 바깥 오차보다 **더 빡세게** (예: `outer_tol * 0.2`)  
- **특이점/NaN 방어**:  
  - `Evaluate()` 실패 시 조기 중단 또는 세분화 재시도  
  - integrand 계산 후 `std::isfinite` 체크 권장

---

## 7) 구현 메모 (코드와 대응)

- `Trap(e)zoid(a,b,j,st)` : 위 “재귀형 사다리꼴”과 동일  
- `QRomberg()` : `K_ROMBERG`개 최근 값을 **1‑based** 배열로 보간 → `h[j-K+1..j]`, `s[j-K+1..j]` 전달  
- `PolynomialInterpolation()` : Neville 보간 기반 외삽 (입력 배열 1‑based 주의)  
- `QSimpsons()` : $\(S_j=\frac{4T_j-T_{j-1}}{3}\)$ 재사용 + 수렴 가드  
- `IsoCurveLengthIntegrator::Evaluate()` : $\(\|S_u\times S_v\|\)$ integrand

---

## 9) 소스 코드

```cpp
class ON_CLASS ON_IntegFuncEvalObj
{
public:
  ON_IntegFuncEvalObj() {}
  virtual ~ON_IntegFuncEvalObj() {}
  virtual bool Evaluate(double dT, double& rdResult) const;
};

inline bool ON_IntegFuncEvalObj::Evaluate(double, double&) const
{
  ON_ASSERT(0); return false;
}

class ON_CLASS ON_IntegratorEx
{
protected:
  const ON_IntegFuncEvalObj&  m_crFunctionEvaluator{};
  double                      m_dS{ 0 };
  double                      m_dA{ 0 };
  double                      m_dB{ 0 };
  IntegratorAlgorithmType     m_eIntegAlgorithm{ IntegratorAlgorithmType::IA_SIMPSONS };
  double                      m_dDesiredAccuracy{ON_TOL6};

public:
  ON_IntegratorEx(const ON_IntegFuncEvalObj& crFunctionEvaluator);
  ~ON_IntegratorEx() {}

  bool IntegrateIt(double dA, double dB,
    IntegratorAlgorithmType  eIntegAlgorithm,
    double dDesiredAccuracy,
    double& rdResult);

  bool PolynomialInterpolation(double* dXA, double* dYA, int lN,
    double dX, double& rdY, double& rdDY);

protected:
  bool Trapezoid(double dA, double dB, int lN, double& rdResult);
  bool QRomberg(double& rdResult);
  bool QSimpsons(double& rdResult);

};

inline ON_IntegratorEx::ON_IntegratorEx(
  const ON_IntegFuncEvalObj& crFunctionEvaluator)
  : m_crFunctionEvaluator(crFunctionEvaluator)
  , m_dS(0.0)
{
}


ON_DECL
double ON_CalcCurveLengthByTangent(
  ON_Curve& curve,
  double dDesiredAccuracy
);


double ON_CalcSurfaceArea(
  const ON_Surface& surface,
  const ON_Interval uDomain,
  const ON_Interval vDomain,
  double dPatchAccuracy
);


#define REFINE_MAX 20         
#define REFINE_MAXP REFINE_MAX+1

#define K_ROMBERG 4

bool ON_IntegratorEx::IntegrateIt(
  double dA, 
  double dB,
  IntegratorAlgorithmType  eIntegAlgorithm,
  double dDesiredAccuracy,
  double& rdResult)
{
  bool ok = false;
  m_dS = 0.0;
  m_dA = dA;
  m_dB = dB;
  m_eIntegAlgorithm = eIntegAlgorithm;
  m_dDesiredAccuracy = dDesiredAccuracy;
  if (eIntegAlgorithm == IntegratorAlgorithmType::IA_SIMPSONS) {
    ok = QSimpsons(rdResult);
  }
  else if (eIntegAlgorithm == IntegratorAlgorithmType::IA_ROMBERG) {
    ok = QRomberg(rdResult);
  }
  else 
  { 
    return false; 
  }

  return ok;
}

#define MAX_INTER 1000

bool ON_IntegratorEx::PolynomialInterpolation(
  double* dXA, 
  double* dYA,
  int iN, 
  double dX,
  double& rdY, 
  double& rdDY)
{
  ON_ASSERT(iN <= MAX_INTER);
  int ns = 1;

  double dif = std::fabs(dX - dXA[0]);

  double c[MAX_INTER];
  double d[MAX_INTER];

  // First find the index ns of the closest table entry
  for (int i = 1; i <= iN; i++) {
    double dift;
    if ((dift = std::fabs(dX - dXA[i])) < dif) {
      ns = i;
      dif = dift;
    }
    c[i] = dYA[i];  // Initialize the tables
    d[i] = dYA[i];
  }

  rdY = dYA[ns--];  // Initial approximation
  // Update things
  for (int m = 1; m < iN; m++) {
    for (int ii = 1; ii <= iN - m; ii++) {
      double ho = dXA[ii] - dX;
      double hp = dXA[ii + m] - dX;
      double w = c[ii + 1] - d[ii];
      double den;
      // This error can occur if two input XA's are identical
      if ((den = ho - hp) == 0.0) return false;
      den = w / den;
      d[ii] = hp * den;
      c[ii] = ho * den;
    }
    if (2 * ns < (iN - m)) {
      rdDY = c[ns + 1];
    }
    else {
      rdDY = d[ns--];
    }
    rdY += rdDY;
  }

  return true;
}

bool ON_IntegratorEx::Trapezoid(double dA, double dB, int iN,
  double& rdResult)
{
  if (iN == 1) { // Compute simple average of end points for initial value
    double dResultA, dResultB;
    m_crFunctionEvaluator.Evaluate(dA, dResultA);
    m_crFunctionEvaluator.Evaluate(dB, dResultB);
    m_dS = 0.5 * (dB - dA) * (dResultA + dResultB);
    rdResult = m_dS;
    return true;
  }

  int it = 1;
  for (int j = 1; j < iN - 1; j++) it <<= 1;
  double tnm = it;
  double del = (dB - dA) / tnm; // Spacing of points to be added
  double x = dA + 0.5 * del;
  double sum = 0.0;
  for (int jj = 1; jj <= it; jj++, x += del) {
    double dResult;
    m_crFunctionEvaluator.Evaluate(x, dResult);
    sum += dResult;
  }
  m_dS = 0.5 * (m_dS + (dB - dA) * sum / tnm);
  rdResult = m_dS;
  return true;
}


bool ON_IntegratorEx::QRomberg(double& rdResult)
{

  double s[REFINE_MAXP + 2], h[REFINE_MAXP + 2];
  double a = m_dA;
  double b = m_dB;

  h[1] = 1.0;
  for (int j = 1; j <= REFINE_MAX; j++) {
    Trapezoid(a, b, j, s[j]);
    if (j >= K_ROMBERG) {
      double ss, dss;
      PolynomialInterpolation(&h[j - K_ROMBERG + 1], &s[j - K_ROMBERG + 1], K_ROMBERG, 0.0, ss, dss);
      if (std::fabs(dss) < m_dDesiredAccuracy * std::fabs(ss)) {
        rdResult = ss;
        return true;
      }
    }
    s[j + 1] = s[j];
    h[j + 1] = 0.25 * h[j];
  }
  return false;
}

bool ON_IntegratorEx::QSimpsons(double& rdResult)
{
  double a = m_dA;
  double b = m_dB;
  double ost = -1.0e30;
  double os = -1.0e30;
  for (int j = 1; j <= REFINE_MAX; j++) {
    double st;
    Trapezoid(a, b, j, st);
    double s = (j == 1) ? st : (4.0 * st - ost) / 3.0;
    if (j > 2)
    {
      if (std::fabs(s - os) < m_dDesiredAccuracy * std::fabs(os)) {
        rdResult = s;
        return true;
      }
    }
    if (s == 0.0 && os == 0.0 && j > 6) {
      rdResult = s;
      return true;
    }
    os = s;
    ost = st;
  }
  return false;
}

class ON_CurveLength : public ON_IntegFuncEvalObj
{
private:
  const ON_Curve& m_crCurve;
public:
  ON_CurveLength(const ON_Curve& crCurve)
    : m_crCurve(crCurve) {}
  bool Evaluate(double dT, double& rdResult) const override;
};

bool ON_CurveLength::Evaluate(double dT, double& rdResult) const
{
  TArrayd evalData(2 * 3);
  m_crCurve.Evaluate(dT, 1, 3, evalData.GetData());
  rdResult = ON_3dVector(&evalData[3]).Length();
  return true;
}

double ON_CalcCurveLengthByTangent(
  ON_Curve& curve, 
  double dDesiredAccuracy
)
{
  double totalLength = 0;
  if (!curve.IsValid()) return totalLength;

  double dSpanAccuracy = dDesiredAccuracy / curve.SpanCount();

  ON_Interval sSegIvl = curve.Domain();
  ON_CurveLength sEval(curve);
  ON_IntegratorEx sIntegrator(sEval);

  double dRomberg;
  sIntegrator.IntegrateIt(sSegIvl.Min(), sSegIvl.Max(), 
    IntegratorAlgorithmType::IA_ROMBERG,
    dSpanAccuracy, dRomberg);
  totalLength += dRomberg;

  return totalLength;
}

class CalcAreaIntegrator : public ON_IntegFuncEvalObj
{
private:
  const ON_Surface&   m_crSurface;
  const ON_Interval&  m_crUDomain;
  const ON_Interval&  m_crVDomain;
  double              m_dIsoCurveAccuracy;
public:
  CalcAreaIntegrator(
    const ON_Surface& crSurface,
    const ON_Interval& crUDomain, 
    const ON_Interval& crVDomain, 
    double dIsoCurveAccuracy)
    : m_crSurface(crSurface), m_crUDomain(crUDomain), m_crVDomain(crVDomain),
    m_dIsoCurveAccuracy(dIsoCurveAccuracy) {}
  bool Evaluate(double dT, double& rdResult) const override;
};


class IsoCurveLengthIntegrator : public ON_IntegFuncEvalObj
{
private:
  const ON_Surface& m_crSurface;
  double m_dUIsoValue;
public:
  IsoCurveLengthIntegrator(const ON_Surface& crSurface, double dUIsoValue)
    : m_crSurface(crSurface), m_dUIsoValue(dUIsoValue) {
  }
  virtual bool Evaluate(double dT, double& rdResult) const;
};


bool IsoCurveLengthIntegrator::Evaluate(double dT, double& rdResult) const
{
  TArrayd taVal(3 * 3);
  m_crSurface.Evaluate(m_dUIsoValue, dT, 1, 3, taVal.GetData());
  ON_3dVector Su = ON_3dVector(&taVal[3]);
  ON_3dVector Sv = ON_3dVector(&taVal[6]);

  double dE = Su.Dot(Su);
  double dF = Su.Dot(Sv);
  double dG = Sv.Dot(Sv);

  rdResult = std::sqrt((std::max)(0.0, dE * dG - dF * dF));
  return true;
}

bool CalcAreaIntegrator::Evaluate(double dT, double& rdResult) const
{
  IsoCurveLengthIntegrator sEval(m_crSurface, dT);
  ON_IntegratorEx sIntegrator(sEval);

  double dRomberg;
  sIntegrator.IntegrateIt(m_crVDomain.Min(),
    m_crVDomain.Max(), IntegratorAlgorithmType::IA_ROMBERG,
    m_dIsoCurveAccuracy, dRomberg);

  rdResult = dRomberg;
  return true;
}

double ON_CalcSurfaceArea(
  const ON_Surface& surface,
  const ON_Interval uDomain,
  const ON_Interval vDomain,
  double dPatchAccuracy)
{
  double totalArea = 0.0;
  CalcAreaIntegrator sAreaEval(surface, uDomain, vDomain, dPatchAccuracy * 0.2);
  ON_IntegratorEx sIntegrator(sAreaEval);
  double dRomberg = 0.0;
  sIntegrator.IntegrateIt(uDomain.Min(),
    uDomain.Max(), IntegratorAlgorithmType::IA_ROMBERG,
    dPatchAccuracy, dRomberg);
  totalArea += dRomberg;

  return totalArea;
}
```


## 8) 짧은 사용 예 (README 데모)

```cpp
// Curve length
double L = ON_CalcCurveLengthByTangent(curve, 1e-8);
printf("length ≈ %.12g\n", L);

// Surface area (full domain)
double A = ON_CalcSurfaceArea(srf, srf.Domain(0), srf.Domain(1), 1e-6);
printf("area ≈ %.12g\n", A);
```
---








