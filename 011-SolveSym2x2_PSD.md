# 📐 SolveSym2x2_PSD(Positive Semi-Definite) 이론 정리

## 1. 약어 의미
- **PSD** = **Positive Semi-Definite**
  - **Positive (양의)**  
  - **Semi (반-, 완전은 아니고 부분적)**  
  - **Definite (정정, definite quadratic form)**  

👉 따라서 PSD는 **양의 준정정 행렬 (양의 반정 행렬)**을 뜻함.

---

## 2. 정의

어떤 대칭 행렬 $\(A\)$ 가 있을 때, 모든 벡터 $\(x \in \mathbb{R}^n\)$ 에 대해

$$
x^T A x \geq 0
$$

가 성립하면 \(A\)는 **Positive Semi-Definite (PSD)** 라고 한다.

---

## 3. PD와 PSD 비교

- **Positive Definite (PD, 양의 정정)**

$$
  x^T A x > 0 \quad \text{(모든 } x \neq 0\text{)}
$$

- **Positive Semi-Definite (PSD, 양의 준정정)** 

$$
  x^T A x \geq 0 \quad \text{(0이 될 수도 있음)}
$$

---

## 4. 요약
- **PD** → 모든 방향에서 양수 (엄밀히 볼록)  
- **PSD** → 양수 또는 0 (볼록하지만 평평한 방향 허용)

## 1. 문제 정의
우리는 종종 **대칭 2×2 행렬**을 다루게 된다:

$$
A =
\begin{bmatrix}
a & b \\
b & c
\end{bmatrix}
$$

이 행렬은 **대칭(symmetric)**이고, 어떤 상황에서는 **반정(positive semidefinite, PSD)** 특성을 가진다.  
즉, 모든 벡터 $\(x \in \mathbb{R}^2\)$ 에 대해

$$
x^T A x \geq 0
$$

이 성질은 최적화, 곡면 근사, Newton/LM 보정에서 매우 중요하다.

---

## 2. 일반적인 해법
대칭 2×2 행렬의 해석은 간단하다.

- **행렬식**:
 
$$
\det(A) = ac - b^2
$$

- **고유값**:
  
$\lambda_{1,2} = \frac{a+c}{2} \;\pm\; \sqrt{\left(\frac{a-c}{2}\right)^2 + b^2}$

- **PSD 조건**:
  1. $\(a \geq 0\)$  
  2. $\(c \geq 0\)$  
  3. $\(\det(A) \geq 0\)$

세 조건이 모두 만족되면, $\(A\)$ 는 PSD다.

---

## 3. 선형 시스템 해법

우리는 보통 다음 시스템을 풀고 싶다:

<p align="center">
  <img src="https://math.vercel.app/?from=A%5Cbegin%7Bbmatrix%7Dx%5C%5Cy%5Cend%7Bbmatrix%7D%3D%5Cbegin%7Bbmatrix%7Df%5C%5Cg%5Cend%7Bbmatrix%7D" />
</p>




즉,

$$
a x + b y = f, \quad
b x + c y = g
$$

### 3.1 일반 역행렬

만약 $\(\det(A) \neq 0\)$ :

<p align="center">
  <img src="https://math.vercel.app/?from=%5Cbegin%7Bbmatrix%7Dx%5C%5Cy%5Cend%7Bbmatrix%7D%3D%5Cfrac%7B1%7D%7Bac-b%5E2%7D%5Cbegin%7Bbmatrix%7Dc%26-b%5C%5C-b%26a%5Cend%7Bbmatrix%7D%5Cbegin%7Bbmatrix%7Df%5C%5Cg%5Cend%7Bbmatrix%7D" />
</p>

---

## 4. PSD 상황에서의 주의점
PSD 행렬은 **반정**이므로, **특이(singular)** 한 경우가 존재할 수 있다.  
즉, $\(\det(A) = 0\)$ 일 수 있으며, 이 경우 고유값 중 하나가 0이 된다.

### 4.1 특이한 경우
- 만약 $\(a=c=0, b=0\) → \(A=0\)$ .  
  이때는 $\(f=g=0\)$ 이면 해가 무수히 많고, 아니면 해가 없다.

- 만약 $\(ac=b^2\)$ → 순위 1(rank-1) 행렬.  
  이때 해는 **직선 공간** 상에 존재하며, 보통 최소제곱 혹은 Moore–Penrose 의사역행렬을 써야 한다.

---

## 5. 의사역행렬(Pseudoinverse) 접근
의사역행렬 $\(A^+\)$ 는 다음 성질을 만족:

$$
x^* = A^+ \, f
$$

이 $\(x^*\)$ 는 최소제곱 해를 의미하며,  
**PSD 특이 상황에서도 안정적인 해**를 준다.

2×2 경우에는 다음과 같이 구할 수 있다:

- 고유분해 $\(A = Q \Lambda Q^T\)$  
- $\(\Lambda^+\)$ 는 0이 아닌 고유값을 역수로 치환  
- $\(A^+ = Q \Lambda^+ Q^T\)$

---

## 6. SolveSym2x2_PSD 절차 요약
1. $\(a,c,b\)$ 를 읽는다.  
2. 판정:
   - $\(\det > \varepsilon\)$ → 정상 역행렬 사용  
   - $\(\det \approx 0\)$ → 의사역행렬로 최소제곱 해 반환  
   - $\(A=0\) → \(f=g=0\)$ 이면 자유도 무한, 아니면 불가능  
3. 해를 반환하면서, 필요하면 `ok` flag로 유효성 표시.

---

## 7. 활용 맥락
- **LM/뉴턴 최적화**: $\(J^T J\)$ 가 PSD → `SolveSym2x2_PSD`로 안정적 해 구함.  
- **곡면 근사**: 2×2 Hessian 행렬 다루는 경우.  
- **수치 안전성**: rank-1 근처에서 바로 explode 하지 않고 graceful degradation.

---

## 8. 결론
`SolveSym2x2_PSD`는 작은 크기(2×2)의 대칭 PSD 행렬을 안정적으로 푸는 유틸리티다.  
- **정규 경우**: 역행렬로 바로 풀이  
- **특이 경우**: 의사역행렬로 최소제곱 해 반환  
- **0 행렬**: 특수 케이스 처리  

이 접근은 **루프 안에서 반복적으로** 호출되더라도 수치적으로 안전하게 수렴을 도와준다.


---


# 📐 Positive Definite (PD) vs Positive Semi-Definite (PSD)

## 1. 정의

- **Positive Definite (PD, 양의 정정)**  
  대칭 행렬 $\(A\)$ 가 모든 **0이 아닌** 벡터 $\(x\)$ 에 대해

$$
  x^T A x > 0
$$

  를 만족할 때, $\(A\)$ 는 PD.

- **Positive Semi-Definite (PSD, 양의 준정정)**  
  대칭 행렬 $\(A\)$ 가 모든 벡터 $\(x\)$ 에 대해

$$
  x^T A x \geq 0
$$

  를 만족할 때, $\(A\)$ 는 PSD.

👉 즉, **PD는 항상 양수**, **PSD는 양수 또는 0**.

---

## 2. 고유값 관점

- **PD**: 모든 고유값이 **엄밀히 양수**  
- **PSD**: 모든 고유값이 **0 이상**

---

## 3. 예제

- **PD 예제**

$$
  A =
  \begin{bmatrix}
  2 & 0 \\
  0 & 3
  \end{bmatrix}
$$

  - 고유값 = 2, 3  
  - 모두 > 0 → **PD**

- **PSD 예제**

$$
  B =
  \begin{bmatrix}
  1 & 0 \\
  0 & 0
  \end{bmatrix}
$$

  - 고유값 = 1, 0  
  - 모두 ≥ 0 → **PSD (PD 아님)**

---

## 4. 기하학적 의미

- **PD**  
  - $\(x^T A x\)$ 가 모든 방향에서 양수  
  - **타원체(ellipsoid)** 정의  
  - Newton 최적화에서 **엄밀히 볼록(convex)**인 함수의 Hessian에 해당

- **PSD**  
  - $\(x^T A x\)$ 가 음수는 없지만, 특정 방향에서는 0  
  - **납작한 타원통(flat direction)** 정의  
  - Newton 최적화에서 **볼록(convex)하지만 평평한 축 존재**

---

## 5. 차이 요약

| 구분 | PD (Positive Definite) | PSD (Positive Semi-Definite) |
|------|-------------------------|-------------------------------|
| 부등식 | $\(x^T A x > 0\)$ (모든 $\(x \neq 0\)$ ) | $\(x^T A x \geq 0\)$ (모든 $\(x\)$ ) |
| 고유값 | 전부 > 0 | 전부 ≥ 0 |
| 행렬식 | > 0 (모든 주 소행렬) | ≥ 0 (모든 주 소행렬) |
| 예시 | diag(2,3) | diag(1,0) |
| 의미 | 완전히 양의 방향만 존재 | 어떤 축은 0 (평평함) |

---

## 6. 한 줄 요약
- **PD** → 완전히 “볼록”한 행렬 (모든 방향에서 양수)  

- **PSD** → “볼록하긴 하지만 일부 평평한 방향”이 존재할 수 있음


















