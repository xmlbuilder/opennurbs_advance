# 🎯 핵심 개념: Knot Vector가 수식에 미치는 영향
## 1. Knot Vector란?
- B-spline이나 NURBS에서 곡선을 정의하는 매개변수 분할 지점입니다.
- 예: [0, 0, 0, 1, 2, 3, 4, 4, 4] 같은 벡터는 degree 2 곡선의 시작과 끝을 고정하는 역할을 합니다.

## 📐 Basis Function 수식이 달라지는 이유
### ✅ Basis Function은 Knot에 따라 정의됨
- Basis function $N_{i,p}(t)$ 는 재귀적으로 정의되며, 다음과 같은 수식을 따릅니다:


$N_{i,0}(t)$ = 
  if $u_i$ <= t < $u_{i+1}$ then 1
  else 0

$N_{i,p}(t)$ = 
  $(t - u_i)$ / $(u_{i+p} - u_i)$ * $N_{i,p-1}(t)$ +
  ($u_{i+p+1}$ - t) / ($u_{i+p+1}$ - $u_{i+1}$) * $N_{i+1,p-1}(t)$



- 여기서 u_i는 Knot vector의 i번째 값입니다.
- 즉, Knot 값이 바뀌면 Basis function의 정의 구간과 형태가 완전히 달라집니다.

## 🔍 FindSpan 수식이 달라지는 이유
### ✅ FindSpan은 Knot vector에서 t가 속한 구간을 찾는 함수
- 일반적인 구현:

```cpp
FindSpan(n, p, t, U):
  if t == U[n+1]: return n
  binary search to find i such that U[i] <= t < U[i+1]
  return i
```

- 여기서 U는 Knot vector, n은 control point 개수 - 1, p는 degree
- Knot vector가 uniform인지, non-uniform인지, clamped인지에 따라
FindSpan의 결과가 달라지고, 그에 따라 Basis function의 활성 영역도 달라집니다.

## 🧠 요약
# 📐 Knot Vector가 Basis & FindSpan 수식에 미치는 영향

| 항목           | 설명                                                                 |
|----------------|----------------------------------------------------------------------|
| Knot Vector    | 곡선의 정의 구간을 결정하는 매개변수 분할 지점                        |
| Basis Function | Knot 값에 따라 정의 구간과 형태가 달라짐                              |
| 수식 구조      | 재귀적 정의: N_{i,p}(t)는 u_i, u_{i+1}, ..., u_{i+p+1}에 따라 결정됨   |
| 영향           | Knot이 바뀌면 Basis의 활성 영역, 연속성, 형태가 모두 달라짐            |


| 항목           | 설명                                                                 |
|----------------|----------------------------------------------------------------------|
| FindSpan       | Knot vector에서 t가 속한 구간(span)을 찾는 함수                       |
| 수식 구조      | 이진 탐색 또는 선형 탐색으로 U[i] <= t < U[i+1] 조건 만족하는 i 반환   |
| 영향           | Knot 분포에 따라 span 위치가 달라지고, Basis 활성 인덱스도 달라짐      |


| 항목           | 결과 요약                                                            |
|----------------|----------------------------------------------------------------------|
| Knot 변경 시   | Basis function 수식과 FindSpan 결과 모두 달라짐                       |
| 실질적 영향    | 곡선의 형태, 연속성, 평가 위치가 바뀌며 렌더링 결과에 직접적인 영향    |



## ✨ 비유로 설명하자면…
Knot vector는 지도고,
FindSpan은 내가 지금 어느 동네에 있는지 찾는 GPS,
Basis function은 그 동네에서 어떤 건물들이 보이는지 결정하는 시야각이에요.
Knot이 바뀌면 지도가 바뀌고, 그에 따라 GPS도 다르게 작동하고, 보이는 풍경도 달라지는 거죠.
---
