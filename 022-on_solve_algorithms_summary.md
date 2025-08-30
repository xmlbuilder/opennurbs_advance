
# on_solve 알고리즘/수식 정리

지금까지 단계별로 적용된 수식과 알고리즘을 요약했습니다.

---

## 1. Curve–Curve (CC) 교차

- **목표**: 두 곡선 `Ca(t_a)`, `Cb(t_b)`의 교차점 찾기.

### 식
- 잔차:  

```math
r = C_a(t_a) - C_b(t_b)
```

- 야코비안(3×2):  

```math
J = [C_a'(t_a), -C_b'(t_b)]
```

### 커플드 뉴턴 보정
- Normal equations (Levenberg 감쇠):  

```math
(J^T J + \lambda I)\,\Delta = -J^T r
```

- 업데이트:  

```math
t_a \leftarrow t_a + \Delta_0, \quad t_b \leftarrow t_b + \Delta_1
```

### 가속
- 곡선 B를 세그먼트 bbox(RTree)로 분할 → 근방 세그먼트만 검사.

---

## 2. Curve–Surface (CS) 교차

- **목표**: 곡선 `C(t)`와 표면 `S(u,v)`의 교차점 찾기.

### 식
- 잔차:  

```math
r = C(t) - S(u,v)
```

- 야코비안(3×3):  

```math
J = [C'(t), -S_u(u,v), -S_v(u,v)]
```

### 커플드 뉴턴 보정
- Normal equations:  

```math
(J^T J + \lambda I)\,\Delta = -J^T r
```

- 업데이트:  

```math
t \leftarrow t + \Delta_t, \quad u \leftarrow u + \Delta_u, \quad v \leftarrow v + \Delta_v
```

### 가속
- 표면을 (Nu×Nv) 셀로 분할, 각 bbox를 RTree에 삽입 → 후보 셀만 검사.

---

## 3. Surface–Surface (SS) 교차 (추적 기반)

- **목표**: 두 표면 `S_A(u,v)`, `S_B(p,q)` 의 교차 곡선 추적.

### 식
- 잔차:  

```math
r = S_A(u,v) - S_B(p,q)
```

- 야코비안(3×4):  

```math
J = [S_{Au}, S_{Av}, -S_{Bp}, -S_{Bq}]
```

### 커플드 뉴턴 보정 (최소제곱)
- Normal equations:  

```math
(J^T J + \lambda I)\,\Delta = -J^T r
```

- 업데이트:  

```math
(u,v,p,q) \leftarrow (u,v,p,q) + \Delta
```

### 추적 알고리즘 (Predictor–Corrector)
1. **Seed 수집**: A 표면을 샘플링 후 B에 투영, 거리 ≤ tol 포인트 채택 → 클러스터링으로 중복 제거.  
2. **탱전트 예측**:  

```math
t = rac{n_A 	imes n_B}{\|n_A 	imes n_B\|}
```

3. **예측**:  

```math
C_{pred} = C + t \cdot h
```

4. **보정**: `C_pred`를 양 표면에 최근접 투영, 중점 반복 보정 → 간극 작으면 h 증가, 크면 h 감소.  
5. **양방향 추적**: dir = ±1 두 방향으로 확장, 시작점 근처 되돌아오면 클로즈 처리.

---

## 4. 보조 기법

### RTree 가속
- CurveAccel: 곡선을 N 세그먼트로 쪼개고 bbox를 RTree에 저장.  
- SurfaceAccel: 표면을 (Nu×Nv) 셀로 나누어 bbox를 RTree에 저장.  
- 쿼리 시 후보 영역만 검사 → 성능 향상.

### 샘플링/브라켓팅
- 곡선/표면을 균등 샘플링 → 거리 함수 최소화로 후보 탐색.  
- 브라켓팅: 샘플 사이에서 sign-change/minima를 찾아 뉴턴 초기값 제공.

### 리샘플/스무딩/간소화
- `ResampleUniform`: 일정 간격으로 폴리라인 재구성.  
- `SmoothLaplacian`: 라플라시안 스무딩, 엔드포인트 고정.  
- `SimplifyByAngle`: 작은 각(5° 이하) 포인트 제거.

### 도메인/경계 처리
- 각 파라미터 업데이트 후 `ClampUV()` 적용.  
- 경계 근접 시 스텝 축소/종료 처리.

---

## 5. 프로젝트 기본값 (on_solve/Defaults.h)

```cpp
struct Parameters {
  double tol = 1e-6;
  double join_tol = 1e-3;
  double step_initial = 0.5;
  double step_min     = 1e-3;
  double step_max     = 5.0;
  int    max_steps    = 4000;

  int sample_curve     = 256;
  int seed_grid        = 48;
  int accel_grid       = 48;
  int accel_curve_segs = 256;

  double resample_len     = 0.25;
  int    smooth_iters     = 2;
  double smooth_alpha     = 0.5;
  double simplify_angle_d = 5.0;
};
```

- `on_solve::Defaults()`로 접근, `on_solve::SetDefaults(P)`로 런타임 변경 가능.

---

## 6. 핵심 포인트

- **수치적 안정성**: Levenberg 감쇠 λ 사용.  
- **성능**: RTree + 샘플링으로 초기 후보 제한.  
- **정확도**: 커플드 뉴턴으로 sub-tol 수준까지 보정.  
- **출력 품질**: 후처리(리샘플/스무딩/간소화)로 안정적 폴리라인 확보.  
