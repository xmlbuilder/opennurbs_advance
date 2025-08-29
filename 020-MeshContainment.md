# Mesh Containment (OpenNURBS) — Math & Practical Notes

> 목적: 여러 개의 **watertight ON\_Mesh** 간에 **포함 / touching / 분리**를 안정적으로 판정하고, 대용량(△ 50만\~100만)에서도 실용 속도를 내기 위한 공식/팁 모음.
>
> GitHub 렌더링 호환: 모든 수식은 **LaTeX 블록**과 **ASCII 대안**을 같이 제공합니다.

---

## 0) 표기

* 메쉬: $A, B$ (경계는 $\partial A$).
* 길이 스케일: $\L = \max(\operatorname{diag}(\text{AABB}_A), \operatorname{diag}(\text{AABB}_B)\$.
* 공차: $\tau_{\text{dist}} = 10^{-7}\,L$, $\varepsilon_V = 10^{-12}\,L^3$ (프로젝트 스케일에 맞게 1e−8\~1e−6 조정).

**ASCII**: `tau_dist = 1e-7 * L`, `epsV = 1e-12 * L^3`.

---

## 1) 표면 간 최소거리 $d_{\min}$

두 메쉬 표면(경계) 사이의 최소거리.

$$
 d_{\min}(A,B)
 = \min_{\mathbf{p}\in \partial A,\; \mathbf{q}\in \partial B} \lVert \mathbf{p} - \mathbf{q} \rVert
$$

**ASCII**: `dmin(A,B) = min_{p in ∂A, q in ∂B} ||p - q||`
의미:

* `≈ 0` → **touching** (점/선/면 접촉)
* `> tau_dist` → **분리(disjoint)** 또는 **포함**이지만 표면 간 **갭** 존재

---

## 2) 체적 & 질량중심 (Divergence theorem; 원점 이동 포함)

수치 안정화를 위해 AABB 중심 $\mathbf{O}$로 평행이동 후 계산.
삼각형 $(\mathbf{v}_0,\mathbf{v}_1,\mathbf{v}_2)$에 대해:

$$
V = \frac{1}{6}\sum_{f}\Big((\mathbf{v}_0\times\mathbf{v}_1)\cdot\mathbf{v}_2\Big),\qquad
\mathbf{C} = \mathbf{O} + \frac{1}{24V}\sum_{f}(\mathbf{v}_0+\mathbf{v}_1+\mathbf{v}_2)\,\Big((\mathbf{v}_0\times\mathbf{v}_1)\cdot\mathbf{v}_2\Big)
$$

**ASCII**:

```
V = (1/6) * Σ_f dot(cross(v0, v1), v2)
C = O + (1/(24*V)) * Σ_f (v0+v1+v2) * dot(cross(v0, v1), v2)
```

> 주의: $\mathbf{C}$는 **삼각형 방향(오리엔테이션) 일관성**이 중요. 체적의 절댓값 `|V|`만 필요하다면 방향 혼잡에 덜 민감하지만, $\mathbf{C}$는 일관된 바깥 방향이 더 정확합니다.

---

## 3) 고체각 합 기반 와인딩 넘버 (선택적; 참조용)

점 $\mathbf{P}$에서 삼각형 $(\mathbf{a},\mathbf{b},\mathbf{c})$의 고체각:

$$
\Omega_{\triangle} = 2\,\arctan\frac{\det[\mathbf{a}-\mathbf{P},\;\mathbf{b}-\mathbf{P},\;\mathbf{c}-\mathbf{P}]}{\lVert\mathbf{a}-\mathbf{P}\rVert\,\lVert\mathbf{b}-\mathbf{P}\rVert\,\lVert\mathbf{c}-\mathbf{P}\rVert + (\mathbf{a}-\mathbf{P})\cdot(\mathbf{b}-\mathbf{P})\,\lVert\mathbf{c}-\mathbf{P}\rVert + \cdots}
$$

총합 $\Omega = \sum_f \Omega_{\triangle}$,
와인딩 $w = |\Omega|/(4\pi)$.
판정: `inside` if $w > 0.5 + 10\,\varepsilon$.

**ASCII**:

```
omega_tri = 2*atan2( det(a-P, b-P, c-P), |a-P||b-P||c-P| + dot(a-P,b-P)|c-P| + dot(b-P,c-P)|a-P| + dot(c-P,a-P)|b-P| )
omega = Σ omega_tri
w = |omega| / (4*pi)
inside if w > 0.5 + 10*eps
```

> 실전에서는 \*\*BVH 레이캐스트(짝홀 규칙)\*\*를 기본 PIS로, 고체각은 교차 검증용으로 보조 사용하는 것을 권장.

---

## 4) 레이캐스트 PIS (짝홀 규칙; BVH 가속)

임의 3방향($+x,+y,+z$)으로 광선을 쏴 교차 개수를 센 뒤 **다수결**.

* 한 방향당 교차 개수: $N = |\{\text{tri} : \text{ray} \cap \text{tri}\}|$
* 판정: `inside` if `N` is odd, `outside` if even.

**ASCII**: `inside(P,B) = majority( parity( raycast(B, P, dx) ), parity(...dy), parity(...dz) )`
꼭짓점/엣지 히트는 half-count 규칙 + 소량 지터로 특이성 회피.

---

## 5) AABB 최소거리(박스-박스)

두 AABB $A,B$의 축별 간격을 $\Delta x,\Delta y,\Delta z$라 하면,

$$
 d^2_{\text{AABB}}(A,B) = \Delta x^2 + \Delta y^2 + \Delta z^2
$$

여기서 $\Delta x = \max(0, B_{\min,x}-A_{\max,x},\; A_{\min,x}-B_{\max,x})$ 등.

**ASCII**:

```
dx = max(0, B.min.x - A.max.x, A.min.x - B.max.x)
(similar for y,z)
d2 = dx*dx + dy*dy + dz*dz
```

---

## 6) 분류 규칙 (겹침 없음 가정)

대표 내부점 $\mathbf{P}_A=\mathbf{C}_A$, $\mathbf{P}_B=\mathbf{C}_B$.

```
AinB = inside(P_A, B)
BinA = inside(P_B, A)
if (AinB && !BinA) → A ⊂ B    (dmin<=tau 이면 "A ⊂ B (touching)")
if (!AinB && BinA) → B ⊂ A
if (!AinB && !BinA) → (dmin<=tau ? touching : disjoint)
if (AinB && BinA) → (|V_A|<|V_B| ? A ⊂ B : |V_B|<|V_A| ? B ⊂ A : duplicate/touching)
```

* 방향 모호 시 \*\*체적(|V|)\*\*로 결정 (작은 쪽이 내부).
* `dmin`은 **터치 vs 갭** 구분 및 로깅용.

---

## 7) 파이프라인 요약

1. 각 메쉬 전처리: AABB/diag, 체적·중심, **MeshBVH**(삼각형 AABB → `ON_RTree`).
2. 부피 내림차순 정렬로 후보 부모 압축.
3. 쌍 후보에 대해: AABB 포함 필터 → `dmin`(BVH 듀얼 탐색) → PIS(레이캐스트 다수결) → 규칙에 따라 라벨.
4. 자신을 포함하는 후보 중 **최소 체적** 부모를 선택 → 포함 트리 완성.

복잡도(평균): BVH 구축 $O(T\log T)$, 레이캐스트/PIS $O(\log T + k)$, 후보쌍 압축으로 전체는 거의 준선형.

---

## 8) 대용량(△50만\~100만) 실전 팁

* **BVH 필수**: 거리·레이캐스트 모두 `ON_RTree` 기반 후보 축소. (이미 샘플 반영)
* **조기종료**: `dmin <= tau_dist` 즉시 반환.
* **병렬화**: 메쉬별 전처리(BVH/체적)와 쌍 분류를 스레드 풀로.
* **공차 스케일링**: 모델 단위/스케일에 맞게 `tau_dist`(1e−8\~1e−6 L) 튜닝.
* **클린업**: degenerate tri(면적≈0) 제거, 중복 정점 weld, self-intersection 제거, (가능하면) face orientation 일관화.
* **메모리**: 삼각형 100만개면 tri AABB \~ 수십 MB. 필요 시 `tris`를 생략하고 face 인덱스만 저장해 원본 메쉬에서 좌표 참조.

---

## 9) 테스트 시나리오 체크리스트

* **포함**: 큰 상자 안 작은 상자 — `contains` + `dmin > 0`.
* **touching**: 면/엣지/점 접촉 — `touching` + `dmin ≈ 0`.
* **disjoint**: 충분한 거리 — `disjoint`, `dmin ≫ tau`.
* **원형 ㄷ-쉘** 내부 오브젝트: `disjoint/touching` (PIS가 외부로 판정).
* **동일/중복 쉘**: `duplicate` 처리 (체적 거의 동일 + dmin≈0).

---

## 10) 참고 구현 포인트 (요약)

* `MinDistanceMeshes_RTree(A,B,tau)` — 작은 쪽 tri 기준, 상대 RTree에 반경 확장 검색, 후보 tri만 tri–tri 거리.
* `PointInside_RayCastBVH(B,P)` — 3방향 레이캐스트, parity 다수결, 꼭짓점/엣지 히트는 half-count + 지터.
* 방향 충돌 시 체적 비교로 결정, 포함 여부와 touching 플래그는 **분리 저장**.

---

## 11) ASCII 전용 공식 모음 (GitHub 수식 깨짐 대비)

```
# Volume & Centroid (shifted by O)
V = (1/6) * Σ_f dot(cross(v0, v1), v2)
C = O + (1/(24*V)) * Σ_f (v0+v1+v2) * dot(cross(v0, v1), v2)

# Solid angle winding (optional)
omega_tri = 2*atan2( det(a-P, b-P, c-P), |a-P||b-P||c-P| + dot(a-P,b-P)|c-P| + dot(b-P,c-P)|a-P| + dot(c-P,a-P)|b-P| )
omega = Σ omega_tri
w = |omega| / (4*pi)  -> inside if w > 0.5 + 10*eps

# AABB min distance
 dx = max(0, B.min.x - A.max.x, A.min.x - B.max.x)
 dy = max(0, B.min.y - A.max.y, A.min.y - B.max.y)
 dz = max(0, B.min.z - A.max.z, A.min.z - B.max.z)
 d2 = dx*dx + dy*dy + dz*dz

# Classification logic
AinB = inside(P_A, B)
BinA = inside(P_B, A)
if (AinB && !BinA) -> A inside B
elif (!AinB && BinA) -> B inside A
elif (!AinB && !BinA) -> dmin<=tau ? touching : disjoint
else -> (|V_A|<|V_B| ? A inside B : |V_B|<|V_A| ? B inside A : duplicate/touching)
```

---

**끝.** 필요하면 이 문서에 씬 BVH/병렬화 설계도 추가해 드릴 수 있어요.
