# ✂️ Chapter 2. Trimmed Surface

## 1. 배경 (Background)

NURBS 곡면은 기본적으로 **직사각형 도메인**을 가집니다:

$$
(u,v) \in [u_{\min},u_{\max}] \times [v_{\min},v_{\max}]
$$

하지만 실제 CAD 모델링에서는 **곡선으로 잘려진 영역**(Trimmed Region)만이 유효한 경우가 많습니다.  
예: 원기둥에서 윗면/아랫면을 잘라낸 영역, 구에서 구멍을 뚫은 영역.

👉 따라서 **Trimmed Surface = Base Surface + Trim Curves**

---

## 2. 기본 정의 (Definition)

- **Base Surface**: $S(u,v)$, NURBS 곡면  
- **Trimming Loops**:
  - **Outer loop**: 면의 외곽을 정의  
  - **Inner loop(s)**: 구멍(hole) 또는 잘려진 영역  

Trimmed Surface는 다음과 같이 정의됩니다:

$$
\Omega = \{ (u,v) \in D \mid (u,v) \text{ lies inside outer loop and outside inner loops} \}
$$

---

## 3. Trim Curve 투영 (Projection of Trimming Curve)

3D에서 정의된 트림 곡선 $C(t)$를 UV 도메인에 투영:

$$
(u(t), v(t)) = \Pi_{UV}( C(t) )
$$

여기서 $\Pi_{UV}$는 `ON_Surface::ClosestPoint` 또는 Newton 기반 `ON_ClosestUV` 알고리즘을 사용.  
(※ Analytic Sphere / Cylinder의 경우는 직접 공식으로 변환 가능)

---

## 4. 자료구조 (Data Structure)

```cpp
class ON_TrimmedSurface {
public:
    ON_NurbsSurface base_surface;
    std::vector<ON_Curve*> outer_loops;
    std::vector<std::vector<ON_Curve*>> inner_loops;

    bool IsInside(const ON_2dPoint& uv) const;
    void AddOuterLoop(const std::vector<ON_Curve*>& loop);
    void AddInnerLoop(const std::vector<ON_Curve*>& loop);
};
```

---

## 5. 점 포함 검사 (Point-In-Trim Test)

- **Step 1**: 점 $Q$를 Base Surface에 투영하여 UV 좌표 $(u,v)$ 획득  
- **Step 2**: Outer Loop 안에 있는지 → Winding number / Ray casting  
- **Step 3**: Inner Loop에 속하지 않는지 검사  

```cpp
bool ON_TrimmedSurface::IsInside(const ON_2dPoint& uv) const {
    if (!PointInLoop(uv, outer_loops)) return false;
    for (auto& hole : inner_loops) {
        if (PointInLoop(uv, hole)) return false;
    }
    return true;
}
```

---

## 6. 테셀레이션 (Tessellation)

Trimmed Surface는 **Base Surface → Adaptive Tessellation** 후  
Trim Loop에 의해 잘라냄.

- Adaptive sampling → Curvature 기반 subdivide  
- Loop clipping → Polygon clipping (Sutherland–Hodgman, CDT)

👉 최종적으로 **삼각형 메쉬 + Loop-boundary 보존**

---

## 7. 구현 샘플 (Implementation Sample)

```cpp
// Trimmed Surface 생성
ON_TrimmedSurface tsrf;
tsrf.base_surface = sphere;

// Outer Loop 추가
std::vector<ON_Curve*> outer;
outer.push_back(circle3D.ProjectToUV(sphere));
tsrf.AddOuterLoop(outer);

// Inside 테스트
ON_2dPoint uv(0.5, 0.5);
bool inside = tsrf.IsInside(uv);
```

---

## 8. 수학적 도식 (Figures)

### 그림 1: Base Surface (직사각형 도메인)

```svg
<svg width="300" height="200" xmlns="http://www.w3.org/2000/svg">
  <rect x="40" y="40" width="200" height="120" fill="none" stroke="black" stroke-width="2"/>
  <text x="120" y="30" font-size="14">UV Domain</text>
</svg>
```

### 그림 2: Trim Loops (외곽 + 구멍)

```svg
<svg width="300" height="200" xmlns="http://www.w3.org/2000/svg">
  <rect x="40" y="40" width="200" height="120" fill="none" stroke="black" stroke-width="2"/>
  <circle cx="140" cy="100" r="50" fill="none" stroke="blue" stroke-width="2"/>
  <circle cx="140" cy="100" r="20" fill="none" stroke="red" stroke-width="2" stroke-dasharray="4"/>
  <text x="190" y="95" font-size="12">Outer Loop</text>
  <text x="165" y="135" font-size="12">Inner Loop</text>
</svg>
```

### 그림 3: 최종 Trimmed 영역

```svg
<svg width="300" height="200" xmlns="http://www.w3.org/2000/svg">
  <rect x="40" y="40" width="200" height="120" fill="none" stroke="black" stroke-width="2"/>
  <circle cx="140" cy="100" r="50" fill="lightblue" stroke="blue" stroke-width="2"/>
  <circle cx="140" cy="100" r="20" fill="white" stroke="red" stroke-width="2" stroke-dasharray="4"/>
  <text x="100" y="170" font-size="12">Trimmed Region</text>
</svg>
```

---

## 9. 활용 (Applications)

- Boolean Operations → 교집합/차집합 경계 곡선 생성 후 Trimmed Surface 생성  
- Surface-Surface Intersection → 트림 곡선 집합 반환  
- Meshing → 유한요소해석/렌더링용 경계 삼각화  

---

## 10. 연습 문제 (Exercises)

1. 주어진 원기둥 면을 위/아래 평면으로 잘라 Trimmed Surface를 구성하라.  
2. Sphere에 원형 구멍을 뚫은 Trimmed Surface를 구현하라.  
3. Adaptive tessellation 후 Trimmed Mesh를 시각화하라.  

---
