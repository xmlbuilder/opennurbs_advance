# ON_HermiteCurve — 수식과 알고리즘 정리 (OpenNURBS 기반)

이 문서는 `ON_HermiteCurve`의 수학적 정의, 평가식, Bezier/NURBS 변환, 바운딩 박스, 평면 오프셋 생성 로직을 소스 코드에 기반해 정리한 것입니다.  
파라미터 범위는 자연구간 $[0,1]$을 가정합니다.

---

## 1) 정의와 파라미터

입력 데이터 (Hermite 조건):

- 시작점/끝점: $P_1, P_2 \in \mathbb{R}^d$  (여기서 $d = 2$ 또는 $3$)
- 시작/끝 접선(속도): $D_1, D_2 \in \mathbb{R}^d$
- 자연 구간: $u \in [0,1]$

큐빅 Hermite 곡선은 **power-basis** 형태로 다음과 같이 표현됩니다.

$$
H(u) = A + uB + u^2 C + u^3 D
$$

여기서

$$
A = P_1, \quad 
B = D_1, \quad 
C = -3P_1 - 2D_1 + 3P_2 - D_2, \quad 
D = 2P_1 + D_1 - 2P_2 + D_2
$$

> $C, D$는 생성자에서 미리 계산/보관합니다.

### 도함수

$$
\begin{aligned}
H'(u) &= B + 2uC + 3u^2D \\
H''(u) &= 2C + 6uD \\
H'''(u) &= 6D
\end{aligned}
$$

---

## 2) Bezier 제어점과의 등가 표현

동일 곡선을 **3차 Bezier** 곡선으로 보면, 제어점 $B_0,B_1,B_2,B_3$는 다음과 같습니다.

$$
\begin{aligned}
B_0 &= P_1 \\
B_1 &= P_1 + \tfrac{1}{3}D_1 \\
B_2 &= P_2 - \tfrac{1}{3}D_2 \\
B_3 &= P_2
\end{aligned}
$$

따라서 Hermite $\leftrightarrow$ Bezier 변환은 직접적이며, 본 구현은 이 제어점 4개로 `ON_BezierCurve`를 구성해 NURBS로도 변환합니다.

---

## 3) 평가(Evaluate) 알고리즘

요청 도함수 차수 $k$에 따라 $H(u), H'(u), H''(u), H'''(u)$를 위 식으로 계산합니다.  
3차를 초과한 고차 도함수는 $0$ 벡터입니다.

의사코드:

```
u = parameter
P0 = A + u*(B + u*(C + u*D))        // H(u)
if k >= 1:  P1 = B + u*(2*C + 3*u*D)
if k >= 2:  P2 = 2*C + 6*u*D
if k >= 3:  P3 = 6*D
if k >= 4:  P4.. = 0
```

---

## 4) 바운딩 박스 (Bounding Box)

Hermite를 등가 Bezier로 본 뒤, 제어점 $\{B_0,B_1,B_2,B_3\}$에 (필요 시 변환행렬 적용 후) 박스를 맞춥니다.

$$
B_0=P_1, \quad B_1=P_1+\tfrac{1}{3}D_1, \quad B_2=P_2-\tfrac{1}{3}D_2, \quad B_3=P_2
$$

> Bezier 제어 다각형의 AABB는 곡선을 항상 포함합니다. (일반적으로 더 타이트한 최적 박스는 아니지만, 구현이 단순하고 안전합니다.)

---

## 5) NURBS 형식으로의 변환

위 4개 Bezier 제어점을 사용해 3차 NURBS 곡선을 구성합니다.

- 차수: 3  
- 제어점: $B_0..B_3$ (비가중/균일)  
- 노트벡터 (Bezier 등가): $[0,0,0,0,1,1,1,1]$

구현은 `ON_BezierCurve::GetNurbForm(ON_NurbsCurve&)`를 호출합니다.

---

## 6) 평면 오프셋 Hermite 생성자 (Bezier 기반 오프셋)

특별 생성자

$$
\mathrm{ON\\_HermiteCurve}(d, \text{curve}, n, \delta)
$$


는 기준 Hermite 곡선(`curve`)을 같은 **오프셋 평면 법선** $n$에 대해 거리 $\delta$만큼 오프셋한 곡선을 생성합니다.

### 과정

1. 원 Hermite를 Bezier 제어점 $P_0..P_3$로 변환:

$$
P_0=P_1, \quad P_1=P_1+\tfrac{1}{3}D_1, \quad P_2=P_2-\tfrac{1}{3}D_2, \quad P_3=P_2
$$

2. 제어 다각형 변:

$$
a_0=P_1-P_0, \quad a_1=P_2-P_1, \quad a_2=P_3-P_2, \quad a_3=P_3-P_0
$$

3. 각 끝변에 대해 오프셋 방향(평면 접선) 정의:

$$
a_0^\top = \frac{a_0 \times n}{\lVert a_0 \times n \rVert}, \quad
a_2^\top = \frac{a_2 \times n}{\lVert a_2 \times n \rVert}
$$

4. 세 가지 경우에 따라 오프셋 제어점 $Q_0..Q_3$ 설정:

- **직선형 (투영 상 공선)**:

$$
Q_0=P_0+\delta a_0^\top, \quad
Q_1=P_1+\delta a_0^\top, \quad
Q_2=P_2+\delta a_2^\top, \quad
Q_3=P_3+\delta a_2^\top
$$

- **끝변 평행**:

$$
Q_1=P_1+\delta a_0^\top + \tfrac{8\delta}{3}\tfrac{a_0}{\lVert a_0 \rVert+\lVert a_2 \rVert}, \quad
Q_2=P_2+\delta a_2^\top - \tfrac{8\delta}{3}\tfrac{a_2}{\lVert a_0 \rVert+\lVert a_2 \rVert}
$$

- **일반 케이스**:

$$
V = 2 \tfrac{a_1+a_3}{\lVert a_1+a_3 \rVert} - \tfrac{a_0}{\lVert a_0 \rVert} - \tfrac{a_2}{\lVert a_2 \rVert}
$$

$$
Q_1 = P_1+\delta a_0^\top + \tfrac{4\delta}{3}\tfrac{\langle V,a_2 \rangle}{\langle a_0, a_2^\top \lVert a_2 \rVert \rangle}a_0, \quad
Q_2 = P_2+\delta a_2^\top + \tfrac{4\delta}{3}\tfrac{\langle V,a_0 \rangle}{\langle a_2, a_0^\top \lVert a_0 \rVert \rangle}a_2
$$

5. 오프셋 Bezier 제어점에서 Hermite 파라미터 재구성:

$$
P_1'=Q_0, \quad D_1'=3(Q_1-Q_0), \quad P_2'=Q_3, \quad D_2'=3(Q_3-Q_2)
$$

새 Hermite의 $C',D'$는 1절 공식으로 재계산합니다.

> $a_0,a_2$가 너무 짧거나, $a_0^\top,a_2^\top$가 퇴화(영벡터)하는 경우는 조기 종료합니다.

---

## 7) 유효성 (Validity)

- $(P_1,D_1,P_2,D_2,C,D)$가 유효하고, 차원 $d \in \{2,3\}$이면 유효  
- 오프셋 생성자에서는 분모/외적 크기 등의 수치적 퇴화를 방지하는 검사 포함

---

## 8) 구현 포인트 요약

- Hermite $\leftrightarrow$ Bezier 변환을 적극 활용  
- $C,D$를 미리 저장해 빠른 평가/도함수 계산  
- 평면 오프셋은 제어다각형 기반의 기하식으로 세 가지 경우 분기  
- 자연 구간 $[0,1]$, 필요 시 `ON_Xform`으로 변환 후 AABB 갱신

---

## 9) 참고: 간단한 사용 예

```cpp
// 3D Hermite 만들기
ON_3dPoint  P1(0,0,0), P2(1,1,0);
ON_3dVector D1(1,0,0), D2(0,1,0);
ON_HermiteCurve H(P1,D1,P2,D2,3);

// 평가 (점과 1~3차 도함수)
ON_3dPoint d[4];
H.Evaluate(0.25, 4, true, d);

// Bezier/NURBS 변환
ON_BezierCurve B;
ON_NurbsCurve  N;
H.GetBezierCurve(B);
H.GetNurbForm(N);

// 평면 오프셋 (법선 n, 거리 delta)
ON_3dVector n(0,0,1);
double delta = 0.1;
ON_HermiteCurve Hoff(3, H, n, delta);
```

---
## 10) 소스 코드
```cpp


class ON_CLASS ON_HermiteCurve 
{
private:
  ON_3dPoint        m_vP1;  // Initial point
  ON_3dVector       m_vD1;  // Initial direction
  ON_3dPoint        m_vP2;  // Final point
  ON_3dVector       m_vD2;  // Final direction
  
  ON_3dVector       m_vC;   // C of power basis eqn.
  ON_3dVector       m_vD;   // D of power basis eqn.
  int               m_dim;

public:
  
  ON_HermiteCurve(const ON_3dPoint& crP1,
    const ON_3dVector& crV1,
    const ON_3dPoint& crP2,
    const ON_3dVector& crV2,
    int lDimension = 3);

  ON_HermiteCurve(const ON_2dPoint& crP1,
    const ON_2dVector& crV1,
    const ON_2dPoint& crP2,
    const ON_2dVector& crV2,
    int lDimension = 2);

  ON_HermiteCurve(int lDimension,
    const ON_HermiteCurve& crHermiteToOffset,
    const ON_3dVector& crOffsetPlaneNormal,
    double dOffsetDistance);

  virtual ~ON_HermiteCurve() {}


  bool IsValid() const;

  virtual bool CalculateBoundingBox(ON_BoundingBox& tight_bbox,
    bool bGrowBox,
    const ON_Xform* xform
  ) const;

  void GetBezierPoints(ON_3dPoint& rP1,
    ON_3dPoint& rP2,
    ON_3dPoint& rP3,
    ON_3dPoint& rP4) const;

  virtual ON_Interval GetNaturalInterval() const
  {
    return ON_Interval(0.0, 1.0);
  }

  virtual bool Evaluate(double dParameter,
    int lNumDerivatives,
    bool bFromLeft,
    ON_3dPoint aPointAndDerivatives[]) const;

  virtual bool EvaluatePoint(double dParameter, ON_3dPoint& rPoint) const;

  bool GetBezierCurve(ON_BezierCurve& bezierCurve) const;
  bool GetNurbForm(ON_NurbsCurve& nurbsCurve) const;

};

inline ON_HermiteCurve::ON_HermiteCurve(
  const ON_3dPoint& crP1,
  const ON_3dVector& crD1,
  const ON_3dPoint& crP2,
  const ON_3dVector& crD2,
  int lDimension)
  : m_vP1(crP1), m_vD1(crD1), m_vP2(crP2), m_vD2(crD2), m_dim(lDimension)
{
  // Store C and D of power form
  m_vC = -3.0 * m_vP1 - 2.0 * m_vD1 + 3.0 * m_vP2 - m_vD2;
  m_vD = 2.0 * m_vP1 + m_vD1 - 2.0 * m_vP2 + m_vD2;
}

inline ON_HermiteCurve::ON_HermiteCurve(const ON_2dPoint& crP1,
  const ON_2dVector& crD1,
  const ON_2dPoint& crP2,
  const ON_2dVector& crD2,
  int lDimension)
  : m_vP1(crP1), m_vD1(crD1), m_vP2(crP2), m_vD2(crD2), m_dim(lDimension)
{
  // Store C and D of power form
  m_vC = -3.0 * m_vP1 - 2.0 * m_vD1 + 3.0 * m_vP2 - m_vD2;
  m_vD = 2.0 * m_vP1 + m_vD1 - 2.0 * m_vP2 + m_vD2;
}

inline void ON_HermiteCurve::GetBezierPoints(ON_3dPoint& rP1,
  ON_3dPoint& rP2,
  ON_3dPoint& rP3,
  ON_3dPoint& rP4) const
{
  rP1 = m_vP1;
  rP2 = m_vP1 + m_vD1 / 3.0;
  rP3 = m_vP2 - m_vD2 / 3.0;
  rP4 = m_vP2;
}
```

```cpp
#include "opennurbs.h"

ON_HermiteCurve::ON_HermiteCurve(int lDimension,
  const ON_HermiteCurve& crHermiteToOffset,
  const ON_3dVector& crOffsetPlaneNormal,
  double dOffsetDistance)
  : m_dim(lDimension)
{
  ON_3dVector sN = crOffsetPlaneNormal;
  sN.Unitize();

  // Convert to Bezier form
  ON_3dPoint sP0 = crHermiteToOffset.m_vP1;
  ON_3dPoint sP1 = crHermiteToOffset.m_vP1 + crHermiteToOffset.m_vD1 / 3.0;
  ON_3dPoint sP2 = crHermiteToOffset.m_vP2 - crHermiteToOffset.m_vD2 / 3.0;
  ON_3dPoint sP3 = crHermiteToOffset.m_vP2;

  // Compute a
  ON_3dVector a0 = sP1 - sP0;
  ON_3dVector a1 = sP2 - sP1;
  ON_3dVector a2 = sP3 - sP2;
  ON_3dVector a3 = sP3 - sP0;

  if (a0.LengthSquared() < ON_EPSILON) return;
  if (a2.LengthSquared() < ON_EPSILON) return;

  // Compute a0 Transpose and a2 Transpose
  ON_3dVector a0T = ON_CrossProduct(a0, crOffsetPlaneNormal);
  ON_3dVector a2T = ON_CrossProduct(a2,  crOffsetPlaneNormal);
  if (a0T.LengthSquared() < ON_EPSILON) return;
  if (a2T.LengthSquared() < ON_EPSILON) return;

  a0T.Unitize();
  a2T.Unitize();

  // Test for first case where all points are on same line (relative to offset plane
  // projection.
  double d = dOffsetDistance;
  ON_3dPoint sQ0, sQ1, sQ2, sQ3;
  sQ0 = sP0 + d * a0T;
  sQ3 = sP3 + d * a2T;
  if (fabs(ON_DotProduct(a1, a0T)) < ON_EPSILON && fabs(ON_DotProduct(a2, a0T)) < ON_EPSILON) {
    // Have straight line.
    sQ1 = sP1 + d * a0T;
    sQ2 = sP2 + d * a2T;
  }
  else if (fabs(ON_DotProduct(a2, a0T)) < ON_EPSILON) {
    // Have case where end edges of control polygon are parallel
    sQ1 = sP1 + d * a0T + (8.0 * d / 3.0) * a0 / (a0.Length() + a2.Length());
    sQ2 = sP2 + d * a2T - (8.0 * d / 3.0) * a2 / (a0.Length() + a2.Length());
  }
  else {
    // Have standard Bezier offset case
    // Compute V

    ON_3dVector a1a3 = a1 + a3;
    if (a1a3.LengthSquared() < ON_EPSILON) return;

    ON_3dVector V = 2.0 * (a1 + a3) / a1a3.Length() - a0 / a0.Length() - a2 / a2.Length();
    sQ1 = sP1 + d * a0T + (4.0 * d / 3.0) * ((ON_DotProduct(V, a2)) / (ON_DotProduct(a0, a2T * a2.Length()))) * a0;
    sQ2 = sP2 + d * a2T + (4.0 * d / 3.0) * ((ON_DotProduct(V, a0)) / (ON_DotProduct(a2, a0T * a0.Length()))) * a2;
  }

  m_vP1 = sQ0;
  m_vD1 = 3.0 * (sQ1 - sQ0);
  m_vP2 = sQ3;
  m_vD2 = 3.0 * (sQ3 - sQ2);
  m_vC = -3.0 * m_vP1 - 2.0 * m_vD1 + 3.0 * m_vP2 - m_vD2;
  m_vD = 2.0 * m_vP1 + m_vD1 - 2.0 * m_vP2 + m_vD2;
}

bool ON_HermiteCurve::IsValid() const
{
  if (!m_vP1.IsValid() || !m_vP2.IsValid())
    return false;
  if (!m_vD1.IsValid() || !m_vD2.IsValid())
    return false;
  if (!m_vC.IsValid() || !m_vD.IsValid())
    return false;

  return m_dim == 2 || m_dim == 3;
}

bool ON_HermiteCurve::CalculateBoundingBox(ON_BoundingBox& tight_bbox,
  bool bGrowBox,
  const ON_Xform* xform
) const
{
 
  // Compute Bezier polygon for Hermite and use points to compute bounding
  // boxes.
  ON_3dVector sPnts[4];
  sPnts[0] = m_vP1;
  sPnts[1] = m_vP1 + m_vD1 / 3.0;
  sPnts[2] = m_vP2 - m_vD2 / 3.0;
  sPnts[3] = m_vP2;

  if (bGrowBox && !tight_bbox.IsValid())
  {
    bGrowBox = false;
  }
  
  if (xform != nullptr && !xform->IsIdentity())
  {
    for (int i = 0; i < 4; i++)
    {
      sPnts[i].Transform(*xform);
      tight_bbox.Set(sPnts[i], true);
    }
  }
  return true;
}

bool ON_HermiteCurve::Evaluate(
  double dParameter,
  int lNumDerivatives,
  bool,
  ON_3dPoint aPointAndDerivatives[]) const
{
  // Initialize derivatives requested greater than third
  for (int i = 4; i < lNumDerivatives; i++) {
    aPointAndDerivatives[i].x = 0.0;
    aPointAndDerivatives[i].y = 0.0;
    aPointAndDerivatives[i].z = 0.0;
  }

  const ON_3dVector& A = m_vP1;
  const ON_3dVector& B = m_vD1;
  const ON_3dVector& C = m_vC;
  const ON_3dVector& D = m_vD;

  double u = dParameter;

  // Compute position
  aPointAndDerivatives[0] = A + u * (B + u * (C + u * D));

  if (lNumDerivatives >= 1) {
    aPointAndDerivatives[1] = B + u * (2.0 * C + 3.0 * u * D);
  }
  if (lNumDerivatives >= 2) {
    aPointAndDerivatives[2] = 2.0 * C + 6.0 * u * D;
  }
  if (lNumDerivatives >= 3) {
    aPointAndDerivatives[3] = 6.0 * D;
  }
  return true;
}

bool ON_HermiteCurve::EvaluatePoint(double dParameter, ON_3dPoint& rPoint) const
{
  const ON_3dVector& A = m_vP1;
  const ON_3dVector& B = m_vD1;
  const ON_3dVector& C = m_vC;
  const ON_3dVector& D = m_vD;

  double u = dParameter;
  rPoint = A + u * (B + u * (C + u * D));

  return true;
}

bool ON_HermiteCurve::GetBezierCurve(ON_BezierCurve& bezierCurve) const
{
  if (!IsValid()) return false;

  ON_3dPointArray points(4);
  points.SetCount(4);
  points[0] = m_vP1;
  points[1] = m_vP1 + m_vD1 / 3.0;
  points[2] = m_vP2 - m_vD2 / 3.0;
  points[3] = m_vP2;
  bezierCurve = ON_BezierCurve(points);
  return bezierCurve.IsValid();
}

bool ON_HermiteCurve::GetNurbForm(ON_NurbsCurve& nurbsCurve) const
{
  if (!IsValid()) return false;

  ON_3dPointArray points(4);
  points.SetCount(4);
  points[0] = m_vP1;
  points[1] = m_vP1 + m_vD1 / 3.0;
  points[2] = m_vP2 - m_vD2 / 3.0;
  points[3] = m_vP2;
  ON_BezierCurve bezierCurve = ON_BezierCurve(points);
  bezierCurve.GetNurbForm(nurbsCurve);
  return nurbsCurve.IsValid();
}


```

---

