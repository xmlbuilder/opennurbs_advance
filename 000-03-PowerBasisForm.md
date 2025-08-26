# Power Basis 변환과 활용

## 1. Bezier 곡선과 Basis

Bezier 곡선은 보통 **Bernstein basis** 로 정의됩니다.

$$
B(t) = \sum_{i=0}^n P_i \, b_{i,n}(t)
$$

여기서 $b_{i,n}(t)$ 는 Bernstein 다항식입니다.

---

## 2. Power basis란?

Power basis는 단순한 다항식 형태입니다:

$$
B(t) = a_0 + a_1 t + a_2 t^2 + \cdots + a_n t^n
$$

---

## 3. Power basis를 쓰는 이유

- **미분/적분이 쉬움** → $a_i t^i$ 형태는 도함수 계산이 단순  
- **다항식 연산 통합** → 곡선-곡선, 곡선-표면 교차 계산에서 공통 다항식 연산 사용  
- **계수 비교 용이** → 두 곡선이 같은지, 평행/일치 여부 판단에 유리  
- **Newton iteration과 같은 수치 알고리즘에서 효율적**

즉, Bezier basis를 Power basis로 바꾸면 계산이 더 단순하고, 다양한 기하 알고리즘에 활용할 수 있습니다.

---

## 4. 변환 공식

Bezier → Power 변환은 다음과 같은 일반 행렬 곱으로 표현됩니다:

$$
a = M_{BP} \, P
$$

여기서:

- $a$ : Power basis 계수 벡터  
- $P$ : Bezier control point 벡터  
- $M_{BP}$ : Bernstein → Power basis 변환 행렬

역변환도 존재합니다:

$$
P = M_{PB} \, a
$$

---

## 5. 예시 (3차 Bezier)

3차 Bezier 곡선은:

$$
B(t) = (1-t)^3 P_0 + 3t(1-t)^2 P_1 + 3t^2(1-t) P_2 + t^3 P_3
$$

이를 전개하면 Power basis로 표현할 수 있습니다:

$$
B(t) = (P_0) + (-3P_0 + 3P_1)t + (3P_0 - 6P_1 + 3P_2)t^2 + (-P_0 + 3P_1 - 3P_2 + P_3)t^3
$$

즉, 변환 행렬 $M_{BP}$ 는 다음과 같습니다:

$$
M_{BP} =
\begin{bmatrix}
1 & 0 & 0 & 0 \\
-3 & 3 & 0 & 0 \\
3 & -6 & 3 & 0 \\
-1 & 3 & -3 & 1
\end{bmatrix}
$$

---

## 6. 요약

- Bezier basis는 직관적이지만 계산에는 다소 불편.  
- Power basis는 계산 효율성과 수치적 안정성에서 강점.  
- 변환 행렬을 통해 양방향 변환이 가능하며, CAD/CG 분야에서 교차, 미분, 적분, 곡선 비교에 널리 활용됨.


## 구현 코드
```cpp


int ON_Factorial(int n)
{
  if (n == 0)
    return 1;
  else
    return n * ON_Factorial(n - 1);
}


double ON_Binomial(int n, int i)
{
  return ON_Factorial(n) / (ON_Factorial(i) * ON_Factorial(n - i));
}


ON_Matrix ON_BezierToPowerMatrix(int degree)
{
  const int n = degree;
  ON_Matrix M(n + 1, n + 1);
  M.Zero();
  for (int m = 0; m <= n; ++m)
  {
    for (int i = 0; i <= n; ++i)
    {
      if (m >= i)
      {
        int k = m - i;
        double val = ON_Binomial(n, i) * ON_Binomial(n - i, k) * ((k % 2) ? -1.0 : 1.0);
        M[m][i] = val;
      }
    }
  }
  return M;
}


// --- Power -> Bezier 변환행렬: M_PB = (M_BP)^{-1} ---
ON_Matrix ON_PowerToBezierMatrix(int degree, const ON_Matrix& M_BP)
{
  ON_Matrix Minv = M_BP;
  const bool ok = Minv.Invert(ON_ZERO_TOLERANCE);
  if (!ok)
    ON_ERROR("PowerToBezierMatrix: inversion failed.");
  return Minv;
}

ON_Matrix ON_BezierToPowerMatrixSuperfluous(
  int degree)
{
  ON_Matrix matrix(degree + 1, degree + 1);
  for (int i = 0; i <= degree; i++)
  {
    for (int j = 0; j <= degree; j++)
    {
      matrix[i][j] = 0.0;
    }
  }

  matrix[0][0] = matrix[degree][degree] = 1.0;
  matrix[degree][0] = degree % 2 == 0 ? -1.0 : 1.0;

  double sign = -1.0;
  for (int i = 1; i < degree; i++)
  {
    matrix[i][i] = ON_Binomial(degree, i);
    matrix[i][0] = matrix[degree][degree - 1] = sign * matrix[i][i];
    sign = -sign;
  }

  int k1 = (degree + 1) / 2;
  int pk = degree - 1;
  for (int k = 1; k < k1; k++)
  {
    sign = -1.0;
    for (int j = k + 1; j <= pk; j++)
    {
      matrix[j][k] = matrix[pk][degree - j] = sign * ON_Binomial(degree, k) * ON_Binomial(degree - k, j - k);
      sign = -sign;
    }
    pk = pk - 1;
  }
  return matrix;
}


// --- ON_BezierCurve에서 좌표별 CV 벡터 추출 (비유리 3D만) ---
bool ON_ExtractBezierCVPoints(
  const ON_BezierCurve& bc,
  std::vector<double>& Px,
  std::vector<double>& Py,
  std::vector<double>& Pz)
{
  if (bc.Dimension() != 3 || bc.IsRational()) return false;
  const int order = bc.Order();
  Px.resize(order); Py.resize(order); Pz.resize(order);
  for (int i = 0; i < order; ++i)
  {
    const double* cv = bc.CV(i);
    Px[i] = cv[0]; Py[i] = cv[1]; Pz[i] = cv[2];
  }
  return true;
}


// --- Power 다항식 평가 (3D) ---
ON_3dPoint ON_EvaluatePower3D(
  const std::vector<double>& ax,
  const std::vector<double>& ay,
  const std::vector<double>& az,
  double t)
{
  return ON_3dPoint(ON_HornerAscending(ax, t),
    ON_HornerAscending(ay, t),
    ON_HornerAscending(az, t));
}

// --- 행렬 * 열벡터 ---
std::vector<double> ON_MulMatrix(
  const ON_Matrix& M,
  const std::vector<double>& x)
{
  const int R = M.RowCount(), C = M.ColCount();
  std::vector<double> y(R, 0.0);
  for (int r = 0; r < R; ++r)
  {
    double s = 0.0;
    for (int c = 0; c < C; ++c) s += M[r][c] * x[c];
    y[r] = s;
  }
  return y;
}

```

## 테스트 코드

```cpp
static void TestGeneralDegree(int degree)
{
  const int order = degree + 1;

  // 예제 Bezier (비유리 3D, CV는 임의로 생성)
  ON_BezierCurve bc(3, false, order);
  for (int i = 0; i < order; ++i)
  {
    // 간단한 패턴: x=i, y=(i%2? 2.0:0.0), z=0
    bc.SetCV(i, ON_3dPoint((double)i, (i % 2 ? 2.0 : 0.0), 0.0));
  }

  // CV 추출
  std::vector<double> Px, Py, Pz;
  if (!ON_ExtractBezierCVPoints(bc, Px, Py, Pz))
  {
    std::cout << "Only non-rational 3D Bezier is supported in this demo.\n";
    return;
  }

  // Bezier -> Power (오름차) 계수
  ON_Matrix M_BP = ON_BezierToPowerMatrix(degree);
  std::vector<double> ax = ON_MulMatrix(M_BP, Px);
  std::vector<double> ay = ON_MulMatrix(M_BP, Py);
  std::vector<double> az = ON_MulMatrix(M_BP, Pz);

  // Power -> Bezier (복원) 확인
  ON_Matrix M_PB = ON_PowerToBezierMatrix(degree, M_BP);
  std::vector<double> Px2 = ON_MulMatrix(M_PB, ax);
  std::vector<double> Py2 = ON_MulMatrix(M_PB, ay);
  std::vector<double> Pz2 = ON_MulMatrix(M_PB, az);

  double max_cv_diff = 0.0;
  for (int i = 0; i < order; ++i)
  {
    max_cv_diff = (std::max)(max_cv_diff, std::fabs(Px[i] - Px2[i]));
    max_cv_diff = (std::max)(max_cv_diff, std::fabs(Py[i] - Py2[i]));
    max_cv_diff = (std::max)(max_cv_diff, std::fabs(Pz[i] - Pz2[i]));
  }

  // t 샘플에서 Bezier vs Power 평가 비교
  double max_eval_diff = 0.0;
  for (int k = 0; k <= 10; ++k)
  {
    double t = k / 10.0;
    double buf[9] = {};
    bc.Evaluate(t, 0, 3, buf);
    ON_3dPoint B(buf[0], buf[1], buf[2]);
    ON_3dPoint P = ON_EvaluatePower3D(ax, ay, az, t);
    max_eval_diff = (std::max)(max_eval_diff, (B - P).Length());
  }

  std::cout << "Degree n=" << degree << "\n";
  std::cout << "  max |CV - recon(CV)| = " << max_cv_diff << "\n";
  std::cout << "  max |Bezier - Power| = " << max_eval_diff << "\n";
}

int main() {
  ON::Begin();
  ON_TextLog log;

  SetConsoleOutputCP(CP_UTF8);
  SetConsoleCP(CP_UTF8);
 
  for (int n : {1, 2, 3, 5, 7})
    TestGeneralDegree(n);

  ON::End();
  return 0;
}

```

# 출력 결과

```
Degree n=1
  max |CV - recon(CV)| = 0
  max |Bezier - Power| = 0
Degree n=2
  max |CV - recon(CV)| = 0
  max |Bezier - Power| = 1.11022e-16
Degree n=3
  max |CV - recon(CV)| = 0
  max |Bezier - Power| = 6.28037e-16
Degree n=5
  max |CV - recon(CV)| = 1.42109e-14
  max |Bezier - Power| = 3.9968e-15
Degree n=7
  max |CV - recon(CV)| = 1.13687e-13
  max |Bezier - Power| = 1.95399e-14
```

---
