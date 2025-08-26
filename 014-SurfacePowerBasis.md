# 📘 Bezier Surface → Power Basis 변환

## 1. 배경

- **Bezier Surface (비지어 곡면)** 은 두 변수 \(u, v\)에 대한 곡면으로, 제어점(Control Point)과 Bernstein Basis로 정의된다.
- 그러나 계산 효율성을 위해 Bezier Basis 대신 **Power Basis** (다항식 기저 \(u^i v^j\))로 변환하는 것이 유용하다.
- 변환은 **선형 행렬 변환**과 **텐서곱 구조**를 이용한다.

---

## 2. Bezier Surface 정의

차수 $\(n\)$ (u 방향), $\(m\)$ (v 방향)의 Bezier Surface:

$$
S(u,v) = \sum_{i=0}^{n} \sum_{j=0}^{m} P_{i,j} B_i^n(u) B_j^m(v)
$$

- $\(P_{i,j}\)$ : 제어점 (Control Point)
- $\(B_i^n(u) = \binom{n}{i} (1-u)^{n-i} u^i\)$ : Bernstein 다항식

---

## 3. Power Basis로의 변환

### Bezier → Power (곡선의 경우)

$$
C(t) = \sum_{i=0}^{n} P_i B_i^n(t) = \sum_{k=0}^{n} a_k t^k
$$

→ 변환 행렬:
$$
\mathbf{a} = M_{BP} \cdot \mathbf{P}
$$

### Bezier Surface → Power Basis
두 방향에 대해 독립적으로 변환 후 **텐서곱** 적용:

$$
A = M_{BP}^{(u)} \cdot P \cdot (M_{BP}^{(v)})^T
$$


- $\(P\)$ : 제어점 행렬 ( $\((n+1) \times (m+1)\)$ )  
- $\(M_{BP}^{(u)}\)$ : u 방향 변환 행렬 ( $\((n+1) \times (n+1)\$ ))  
- $\(M_{BP}^{(v)}\)$ : v 방향 변환 행렬 ( $\((m+1) \times (m+1)\$ ))  
- $\(A\)$ : Power Basis 계수 행렬  

---

## 4. 예시 (Bi-Cubic Bezier Surface, 3차)

$$
S(u,v) = \sum_{i=0}^{3}\sum_{j=0}^{3} P_{i,j} B_i^3(u) B_j^3(v)
$$

변환 후:

$$
S(u,v) = \sum_{p=0}^{3} \sum_{q=0}^{3} a_{pq} \, u^p v^q
$$

여기서 계수 \(a_{pq}\)는:

$$
A = M_{BP}^{(u)} \, P \, (M_{BP}^{(v)})^T
$$

---

## 5. 구현 샘플 (OpenNURBS)

```cpp
int degree_u = 3, degree_v = 3;
ON_BezierSurface surf(3, false, degree_u+1, degree_v+1, 4, 4);

// 변환 행렬
ON_Matrix M_u = ON_BezierToPowerMatrix(degree_u);
ON_Matrix M_v = ON_BezierToPowerMatrix(degree_v);

// X좌표 계수 변환
ON_Matrix Px(4,4);
for(int i=0;i<=degree_u;i++){
  for(int j=0;j<=degree_v;j++){
    ON_3dPoint cv; 
    surf.GetCV(i,j,cv);
    Px[i][j]=cv.x;
  }
}
ON_Matrix Ax = M_u * Px * M_v.Transpose();
```
---
