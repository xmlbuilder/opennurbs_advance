# NURBS Geometry Overview (OpenNURBS)

## ✨ What are NURBS?
NURBS (Non-Uniform Rational B-Splines) provide a unified mathematical basis for representing:
- **Analytic shapes** (conics, quadrics, circles, ellipses, etc.)
- **Free-form entities** (car bodies, ship hulls, organic surfaces)

Key features:
- Intuitive and geometric interpretation
- Fast and numerically stable algorithms
- Invariant under common transformations (translation, rotation, projection)
- Generalization of non-rational B-splines and rational/non-rational Bézier

---

## 📐 Components of a NURBS Curve
A NURBS curve is defined by four things:

1. **Degree**  
   - Positive integer (commonly 1, 2, 3, or 5)  
   - Degree 1 → linear, Degree 2 → quadratic, Degree 3 → cubic, Degree 5 → quintic  
   - Order = degree + 1

2. **Control Points**  
   - Minimum (degree + 1) points  
   - Shape controlled by moving points  
   - Each has a **weight** → rational vs non-rational  
   - Circles/ellipses are always rational

3. **Knots**  
   - Knot vector: (N + degree – 1) numbers, N = # control points  
   - Conditions: non-decreasing sequence, multiplicity ≤ degree  
   - Uniform vs non-uniform  
   - Knot multiplicity controls smoothness vs kinks  
   - Adding knots → adds control points (no shape change)  
   - Removing knots → changes shape

4. **Evaluation Rule**  
   - Formula combining degree, control points, knots  
   - Uses **B-spline basis functions**  
   - Input: parameter `u` → Output: 3D point on curve

---

## 🧮 Knot & Control Point Relationship
- **Degree 1 (polyline):** each knot ↔ one control point  
- **Higher degree:** groups of (degree+1) control points ↔ groups of 2×degree knots  

---

## 📚 References
- Böhm, W., Farin, G., Kahmann, J. (1984). *Survey of CAGD methods*  
- De Boor, C. (1978). *A Practical Guide to Splines*  
- Farin, G. (1997). *Curves and Surfaces for CAGD*  

---

# 📖 NURBS란 무엇인가? (KOR)

NURBS(Non-Uniform Rational B-Splines: 비균일 유리 B스플라인)는 단순한 2D 선, 원, 호부터 복잡한 3D 자유 곡면과 솔리드까지 표현 가능한 수학적 방법입니다.

### ✅ 장점
- 산업 표준 교환 포맷 지원 (STEP, IGES 등)
- 수학적으로 잘 정의되어 있고, 전 세계 대학/연구소에서 교육
- 숙련된 프로그래머, 엔지니어, 디자이너들이 활용
- 직선, 원, 타원, 구, 토러스 같은 **정형 기하학**과 자유 곡면 모두 정확하게 표현
- 삼각형 mesh 기반 근사보다 훨씬 적은 데이터로 표현 가능
- 계산 규칙이 효율적이고 컴퓨터 구현이 용이

---

## 📐 NURBS 커브 정의 요소
- **차수 (degree)**: 일반적으로 1, 2, 3, 5  
- **제어점 (control points)**: 최소 degree+1, 가중치 weight 보유  
- **매듭점 (knots)**: (N+degree–1)개 숫자, 균일/비균일, 중복도 → 곡선의 연속성 결정  
- **계산 규칙 (evaluation rule)**: 파라미터 입력 → 점 위치 출력

---

## 📝 용어 요약
| 용어 | 설명 |
|------|------|
| Degree (차수) | 보통 1,2,3,5. 곡선의 곡률 제어 |
| Order (위수) | Degree+1 |
| Control Points | 곡선/곡면 형태를 결정하는 점 집합 |
| Weight | Rational(유리)/Non-rational(비유리) 구분 |
| Knot Vector | 곡선의 분할과 연속성 제어 |
| Basis Functions | B-spline 기저 함수, 곡선 계산 핵심 |

---

# 🔧 Practical Examples with OpenNURBS

### 1. Creating a Simple NURBS Curve
```cpp
#include "opennurbs.h"

// Example: quadratic curve (degree 2)
void CreateNurbsCurve()
{
  const int degree = 2;
  const int dim = 3; // 3D curve
  const bool rational = false;
  const int cv_count = 4;

  ON_NurbsCurve curve(dim, rational, degree+1, cv_count);

  // Define control points
  curve.SetCV(0, ON_3dPoint(0,0,0));
  curve.SetCV(1, ON_3dPoint(2,1,0));
  curve.SetCV(2, ON_3dPoint(4,0,0));
  curve.SetCV(3, ON_3dPoint(6,2,0));

  // Knot vector (must be non-decreasing)
  curve.m_knot[0] = 0.0;
  curve.m_knot[1] = 0.0;
  curve.m_knot[2] = 0.5;
  curve.m_knot[3] = 1.0;
  curve.m_knot[4] = 1.0;

  if(curve.IsValid())
    ON_wString str(L"NURBS curve created successfully");
}
```

---

### 2. Creating a NURBS Surface (bilinear patch)
```cpp
#include "opennurbs.h"

void CreateNurbsSurface()
{
  const int u_degree = 1;
  const int v_degree = 1;
  const int dim = 3;
  const bool rational = false;

  const int u_cv_count = 2;
  const int v_cv_count = 2;

  ON_NurbsSurface srf(dim, rational, u_degree+1, v_degree+1, u_cv_count, v_cv_count);

  // Define control points
  srf.SetCV(0,0, ON_3dPoint(0,0,0));
  srf.SetCV(1,0, ON_3dPoint(10,0,0));
  srf.SetCV(0,1, ON_3dPoint(0,10,0));
  srf.SetCV(1,1, ON_3dPoint(10,10,5));

  // Knot vectors
  srf.m_knot[0][0] = 0.0;
  srf.m_knot[0][1] = 1.0;
  srf.m_knot[1][0] = 0.0;
  srf.m_knot[1][1] = 1.0;

  if(srf.IsValid())
    ON_wString str(L"NURBS surface created successfully");
}
```

---

### 3. Evaluating a NURBS Curve Point
```cpp
ON_3dPoint EvaluateCurveAt(const ON_NurbsCurve& curve, double t)
{
  ON_3dPoint pt;
  ON_3dVector der;
  if(curve.Ev1Der(t, pt, der)) // evaluate point + first derivative
  {
    return pt;
  }
  return ON_3dPoint::UnsetPoint;
}
```

---

# 🚀 Applications in CAD/Robotics
- **CAD Modeling**: precise definition of free-form surfaces (car body, ship hull, etc.)  
- **Robotics**: smooth path planning using NURBS curves (trajectory interpolation)  
- **Engineering Analysis (FEM)**: geometry exchange with IGES/STEP while preserving exact math definitions  
- **Animation/Graphics**: modeling organic characters with accurate curves/surfaces  

---

# 📌 Conclusion
- NURBS are the industry-standard representation for curves and surfaces.  
- OpenNURBS provides a robust API for creating, editing, and evaluating NURBS geometry.  
- By understanding **degree, control points, knots, and evaluation rules**, developers can implement CAD, robotics, or simulation systems with high precision.
