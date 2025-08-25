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


ON_Matrix ON_BezierToPowerMatrixSuperfluous(
  int degree)
{
  ON_Matrix matrix(degree + 1, degree + 1);
  for (int i = 0; i < degree; i++)
  {
    for (int j = i + 1; j <= degree; j++)
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


ON_Matrix ON_PowerToBezierMatrixSuperfluous(
  int degree, 
  const ON_Matrix& matrix)
{
  ON_Matrix inverseMatrix(degree + 1, degree + 1);
  for (int i = 0; i < degree; i++)
  {
    for (int j = i + 1; j <= degree; j++)
    {
      inverseMatrix[i][j] = 0.0;
    }
  }

  for (int i = 0; i <= degree; i++)
  {
    inverseMatrix[i][0] = inverseMatrix[degree][i] = 1.0;
    inverseMatrix[i][i] = 1.0 / (matrix[i][i]);
  }

  int k1 = (degree + 1) / 2;
  int pk = degree - 1;

  for (int k = 1; k < k1; k++)
  {
    for (int j = k + 1; j < pk; j++)
    {
      double d = 0.0;
      for (int i = k; i < j; i++)
      {
        d = d - matrix[j][i] * inverseMatrix[i][k];
      }
      inverseMatrix[j][k] = d / (inverseMatrix[j][j]);
      inverseMatrix[pk][degree - j] = inverseMatrix[j][k];
    }
    pk = pk - 1;
  }
  return inverseMatrix;
}
```

## 구현 코드

```cpp



```

---
