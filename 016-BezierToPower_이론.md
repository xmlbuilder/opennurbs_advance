# Bezier → Power Basis 변환 행렬 비교

## 개요

두 함수는 Bezier basis를 Power basis로 변환하는 행렬을 생성합니다.

- `ON_BezierToPowerMatrixSuperfluous`: De Boor 기반의 정석적인 수식 사용
- `ON_BezierToPowerMatrix`: OpenNURBS에서 사용하는 실용적인 방식

---

## 함수별 수식

### 1. ON_BezierToPowerMatrixSuperfluous (De Boor 기반)

$$
B_i^n(t) = Σ_{j=i}^{n} [ (-1)^{j - i} * C(n, i) * C(n - i, j - i) ] * t^j
$$

- $B_i^n(t)$ : i번째 Bernstein basis 함수
- $C(n, k)$ : 이항계수 (n choose k)
- $(-1)^{j - i}$ : 교차 부호
- $t^j$ : Power basis 항

이 수식은 Bernstein basis를 Power basis로 표현하는 정석적인 방식입니다.

---

### 2. ON_BezierToPowerMatrix (OpenNURBS 방식)


M[m][i] = 
  if m >= i:
    $C(n, i)$ * $C(n - i, m - i)$ * $(-1)^{m - i}$
  else:
    0

- M[m][i]: 변환 행렬의 m행 i열 값
- $C(n, i)$ , $C(n - i, m - i)$ : 이항계수
- $(-1)^{m - i}$ : 부호 처리


이 수식은 Power basis 항을 Bernstein basis로 표현하는 방식입니다.

---

## 수식 방향 차이

| 항목             | Superfluous (De Boor)         | ON_BezierToPowerMatrix (OpenNURBS) |
|------------------|-------------------------------|-------------------------------------|
| 기준             | Bernstein → Power             | Power ← Bernstein                   |
| 수식 방향        | $B_i^n(t)$ → $t^j$                | $t^m$ → $B_i^n(t)$                      |
| 부호 처리        | $(-1)^{j - i}$                  | $(-1)^(m - i)$                        |
| 이항계수 구조    | $C(n, i) * C(n - i, j - i)$     | $C(n, i) * C(n - i, m - i)$           |
| 행렬 의미        | Bernstein basis → Power basis | Power basis → Bernstein basis       |

---

## 민감한 수식 요소

- 이항계수 계산 (C(n, k)): 큰 차수에서는 precision 문제 발생 가능
- 부호 처리 $(-1)^k$ : 부호 하나만 틀려도 결과가 완전히 달라짐
- 인덱스 순서 (j - i, m - i): 인덱스 오류는 전체 행렬을 왜곡시킴
- 대칭성 활용 여부: 잘못 적용 시 항 누락 또는 중복 가능

---

## 결론

- ON_BezierToPowerMatrixSuperfluous는 De Boor의 정석적인 수학적 수식을 따름
- ON_BezierToPowerMatrix는 OpenNURBS에서 실용적으로 쓰이는 방식
- 두 방식 모두 수학적으로 타당하지만, 행렬의 의미와 방향이 다르므로 결과가 다를 수 있음
