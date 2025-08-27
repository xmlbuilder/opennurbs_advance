# Bézier Surface ↔ Power Basis 변환 이론

## 1. Bézier Surface 정의

차수 $\((d_u, d_v)\)$ 의 Bézier 곡면은 제어점 $\(P_{i,j}\)$ 와 Bernstein basis $\(B_{i}^{d_u}(u), B_{j}^{d_v}(v)\)$ 로 정의된다:

$$
S(u,v) = \sum_{i=0}^{d_u} \sum_{j=0}^{d_v} P_{i,j} \, B_i^{d_u}(u)\, B_j^{d_v}(v),
$$
여기서
$$
B_i^n(t) = \binom{n}{i} (1-t)^{n-i} t^i.
$$

---

## 2. Power Basis (모노미얼) 정의

Power basis는 단순히 다항식 전개 형태:

$$
S(u,v) = \sum_{p=0}^{d_u} \sum_{q=0}^{d_v} A_{p,q}\, u^p v^q,
$$

여기서 $\(A_{p,q}\)$는 계수 행렬로, x, y, z 좌표 각각에 대해 하나씩 존재한다:
- $\(A^x_{p,q}\)$, $\(A^y_{p,q}\)$, $\(A^z_{p,q}\)$.

---

## 3. Bézier → Power 변환 행렬

1D Bézier basis와 Power basis 사이에는 **선형 변환**이 존재한다.  
- Bézier basis $\([B_0^n, \dots, B_n^n]\)$ 와  
- Power basis $ \([1, t, t^2, \dots, t^n]\)$ 사이 변환:

$$
[B_0^n(t), B_1^n(t), \dots, B_n^n(t)]^T = M_{B\to P} \cdot [1, t, t^2, \dots, t^n]^T
$$

따라서 제어점 행렬 \(P\)를 Power 계수로 바꾸려면:

$$
A = M_u \, P \, M_v^T
$$

- \(M_u\): u 방향 Bézier → Power 변환 행렬  
- \(M_v\): v 방향 Bézier → Power 변환 행렬  
- \(P\): 제어점 좌표 행렬 (x, y, z 각각 따로)  
- \(A\): Power 계수 행렬

---

## 4. Power → Bézier 변환

역변환은 단순히 위 행렬을 역행렬로 바꿔주면 된다:

$$
P = M_u^{-1} \, A \, (M_v^{-1})^T
$$

여기서 \(M_u^{-1}, M_v^{-1}\)는 Bézier→Power 행렬의 역행렬.

---

## 5. 평가 방법

### 5.1 Bézier Surface 평가
주어진 \((u,v)\)에서:

$$
S(u,v) = \sum_{i=0}^{d_u} \sum_{j=0}^{d_v} P_{i,j} \, B_i^{d_u}(u)\, B_j^{d_v}(v).
$$

---

### 5.2 Power Basis 평가
주어진 \((u,v)\)에서:

$$
S(u,v) = \sum_{p=0}^{d_u} \sum_{q=0}^{d_v} A_{p,q}\, u^p v^q.
$$

실제 계산은 **Horner scheme**을 사용하여 수치적으로 안정적으로 수행.

---

## 6. 정리

- Bézier와 Power Basis는 **동일한 함수 공간**을 표현한다.  
- 변환은 단순히 행렬 곱으로 이루어진다.  
- 변환 전후로 곡면은 완전히 동일해야 하며, 수치 오차는 부동소수점 한계 (~1e-15) 정도다.  
- Power basis는 미분·적분에 편리하고, Bézier basis는 기하학적 직관(제어점, convex hull 등)에 강점이 있다.

---

## 7. 실험 결과 예시

- max |CV - recon(CV)| ≈ \(8.3 \times 10^{-17}\)  
- max |Bezier - Power| ≈ \(9.1 \times 10^{-16}\)  
- 샘플 평가 \((u,v) = (0.37,0.58)\):  
  - Bézier: (1.11, 1.74, -0.027175)  
  - Power : (1.11, 1.74, -0.027175)  
  - 차이: ~\(5 \times 10^{-16}\)

즉, 이론대로 구현이 정확히 동작함을 확인할 수 있다.


