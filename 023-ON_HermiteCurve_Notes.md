# ON_HermiteCurve — 수식과 알고리즘 정리 (OpenNURBS 기반)

> 이 문서는 `ON_HermiteCurve`의 수학적 정의, 평가식, Bezier/NURBS 변환, 바운딩 박스, 평면 오프셋 생성 로직을 소스 코드에 기반해 정리한 것입니다. 파라미터 범위는 자연구간 `[0,1]`을 가정합니다.

---

## 1) 정의와 파라미터

- 입력 데이터 (Hermite 조건):
  - 시작점/끝점: \( \mathbf{P}_1, \mathbf{P}_2 \in \mathbb{R}^d \) (여기서 \(d = 2\) 또는 \(3\))
  - 시작/끝 접선(속도): \( \mathbf{D}_1, \mathbf{D}_2 \in \mathbb{R}^d \)
- 자연 구간: \( u \in [0,1] \)

큐빅 Hermite 곡선은 **power-basis** 형태로 다음과 같이 표현됩니다.

\[
\mathbf{H}(u) \;=\; \mathbf{A} \;+\; u\,\mathbf{B} \;+\; u^2\,\mathbf{C} \;+\; u^3\,\mathbf{D}.
\]

여기서 본 구현은
\[
\mathbf{A}=\mathbf{P}_1,\quad 
\mathbf{B}=\mathbf{D}_1,\quad 
\mathbf{C}=-3\mathbf{P}_1 - 2\mathbf{D}_1 + 3\mathbf{P}_2 - \mathbf{D}_2,\quad 
\mathbf{D}= 2\mathbf{P}_1 + \mathbf{D}_1 - 2\mathbf{P}_2 + \mathbf{D}_2.
\]

> 위 \(\mathbf{C}, \mathbf{D}\)는 생성자에서 미리 계산/보관합니다.

### 도함수
\[
\begin{aligned}
\mathbf{H}'(u) &= \mathbf{B} + 2u\,\mathbf{C} + 3u^2\,\mathbf{D},\\
\mathbf{H}''(u) &= 2\,\mathbf{C} + 6u\,\mathbf{D},\\
\mathbf{H}'''(u) &= 6\,\mathbf{D}.
\end{aligned}
\]

---

## 2) Bezier 제어점과의 등가 표현

동일 곡선을 **3차 Bezier** 곡선으로 보면, 제어점 \(\mathbf{B}_0,\mathbf{B}_1,\mathbf{B}_2,\mathbf{B}_3\)는 다음과 같습니다.

\[
\mathbf{B}_0 = \mathbf{P}_1,\qquad
\mathbf{B}_1 = \mathbf{P}_1 + \frac{1}{3}\mathbf{D}_1,\qquad
\mathbf{B}_2 = \mathbf{P}_2 - \frac{1}{3}\mathbf{D}_2,\qquad
\mathbf{B}_3 = \mathbf{P}_2.
\]

따라서 Hermite \(\leftrightarrow\) Bezier 변환은 직접적이며, 본 구현은 이 제어점 4개로 `ON_BezierCurve`를 구성해 NURBS로도 변환합니다.

---

## 3) 평가(Evaluate) 알고리즘

요청 도함수 차수 \(k\)에 따라 \(\mathbf{H}(u),\mathbf{H}'(u),\mathbf{H}''(u),\mathbf{H}'''(u)\)를 위 식으로 계산합니다. 
3차를 초과한 고차 도함수는 0 벡터입니다.

의사코드:
```
u = parameter
P0 = A + u*(B + u*(C + u*D))        // H(u)
if k >= 1:  P1 = B + u*(2*C + 3*u*D)
if k >= 2:  P2 = 2*C + 6*u*D
if k >= 3:  P3 = 6*D
if k >= 4:  P4.. = 0
```

---

## 4) 바운딩 박스(Bounding Box)

Hermite를 등가 Bezier로 본 뒤, 제어점 \(\mathbf{B}_0..\mathbf{B}_3\)에 **(필요 시 변환행렬 적용 후)** 박스를 맞춥니다.

\[
\mathbf{B}_0=\mathbf{P}_1,\quad 
\mathbf{B}_1=\mathbf{P}_1+\mathbf{D}_1/3,\quad
\mathbf{B}_2=\mathbf{P}_2-\mathbf{D}_2/3,\quad
\mathbf{B}_3=\mathbf{P}_2.
\]

> Bezier 제어 다각형의 AABB는 곡선을 항상 포함합니다. (일반적으로 “tighter” 최적 박스는 아니지만, 구현이 단순하고 안전합니다.)

---

## 5) NURBS 형식으로의 변환

위 4개 Bezier 제어점을 사용해 3차 NURBS 곡선을 구성합니다.
- 차수: 3
- 제어점: \(\mathbf{B}_0..\mathbf{B}_3\) (비가중/균일)
- 노트벡터(Bezier 등가): \([0,0,0,0,1,1,1,1]\)

구현은 `ON_BezierCurve::GetNurbForm(ON_NurbsCurve&)`를 호출합니다.

---

## 6) 평면 오프셋 Hermite 생성자 (Bezier 기반 오프셋)

특별 생성자
\[
\text{ON\_HermiteCurve}(d,\ \text{curve},\ \mathbf{n},\ \delta)
\]
는 기준 Hermite 곡선(`curve`)을 같은 **오프셋 평면 법선** \(\mathbf{n}\)에 대해 거리 \(\delta\)만큼 오프셋한 곡선을 생성합니다. 구현 요지는 다음과 같습니다.

1. 원 Hermite를 **Bezier 제어점** \(\mathbf{P}_0..\mathbf{P}_3\)로 변환.
   \[
   \mathbf{P}_0=\mathbf{P}_1, \ 
   \mathbf{P}_1=\mathbf{P}_1+\mathbf{D}_1/3, \ 
   \mathbf{P}_2=\mathbf{P}_2-\mathbf{D}_2/3, \ 
   \mathbf{P}_3=\mathbf{P}_2.
   \]

2. 제어 다각형 변:
   \[
   \mathbf{a}_0=\mathbf{P}_1-\mathbf{P}_0,\quad
   \mathbf{a}_1=\mathbf{P}_2-\mathbf{P}_1,\quad
   \mathbf{a}_2=\mathbf{P}_3-\mathbf{P}_2,\quad
   \mathbf{a}_3=\mathbf{P}_3-\mathbf{P}_0.
   \]

3. 각 끝변에 대해, 오프셋 방향(평면 접선)으로
   \[
   \mathbf{a}_0^\top = \frac{\mathbf{a}_0 \times \mathbf{n}}{\|\mathbf{a}_0 \times \mathbf{n}\|},\quad
   \mathbf{a}_2^\top = \frac{\mathbf{a}_2 \times \mathbf{n}}{\|\mathbf{a}_2 \times \mathbf{n}\|}.
   \]

4. 세 가지 경우에 따라 오프셋 제어점 \(\mathbf{Q}_0..\mathbf{Q}_3\)를 설정:
   - **직선형(투영 상 공선)**: 
     \(\mathbf{Q}_0=\mathbf{P}_0+\delta \mathbf{a}_0^\top,\ 
       \mathbf{Q}_1=\mathbf{P}_1+\delta \mathbf{a}_0^\top,\ 
       \mathbf{Q}_2=\mathbf{P}_2+\delta \mathbf{a}_2^\top,\ 
       \mathbf{Q}_3=\mathbf{P}_3+\delta \mathbf{a}_2^\top.\)
   - **끝변 평행**: 
     \(\mathbf{Q}_1=\mathbf{P}_1+\delta \mathbf{a}_0^\top + \frac{8\delta}{3}\frac{\mathbf{a}_0}{\|\mathbf{a}_0\|+\|\mathbf{a}_2\|},\ 
       \mathbf{Q}_2=\mathbf{P}_2+\delta \mathbf{a}_2^\top - \frac{8\delta}{3}\frac{\mathbf{a}_2}{\|\mathbf{a}_0\|+\|\mathbf{a}_2\|}.\)
   - **일반 케이스**:
     \[
     \mathbf{V} = 2\,\frac{\mathbf{a}_1+\mathbf{a}_3}{\|\mathbf{a}_1+\mathbf{a}_3\|}
     - \frac{\mathbf{a}_0}{\|\mathbf{a}_0\|}
     - \frac{\mathbf{a}_2}{\|\mathbf{a}_2\|},
     \]
     \[
     \mathbf{Q}_1=\mathbf{P}_1+\delta \mathbf{a}_0^\top + \frac{4\delta}{3}\,\frac{\langle \mathbf{V},\mathbf{a}_2\rangle}{\langle \mathbf{a}_0,\ \mathbf{a}_2^\top \|\mathbf{a}_2\|\rangle}\,\mathbf{a}_0,\quad
     \mathbf{Q}_2=\mathbf{P}_2+\delta \mathbf{a}_2^\top + \frac{4\delta}{3}\,\frac{\langle \mathbf{V},\mathbf{a}_0\rangle}{\langle \mathbf{a}_2,\ \mathbf{a}_0^\top \|\mathbf{a}_0\|\rangle}\,\mathbf{a}_2.
     \]

5. 오프셋 Bezier 제어점 \(\mathbf{Q}_0..\mathbf{Q}_3\)로부터 Hermite 파라미터를 역산:
   \[
   \mathbf{P}_1'=\mathbf{Q}_0,\quad
   \mathbf{D}_1'=3(\mathbf{Q}_1-\mathbf{Q}_0),\quad
   \mathbf{P}_2'=\mathbf{Q}_3,\quad
   \mathbf{D}_2'=3(\mathbf{Q}_3-\mathbf{Q}_2),
   \]
   그리고 새 Hermite의 \(\mathbf{C}',\mathbf{D}'\)는 1절 공식으로 재계산.

> \(\mathbf{a}_0,\mathbf{a}_2\)가 너무 짧거나, \(\mathbf{a}_0^\top,\mathbf{a}_2^\top\)가 퇴화(영벡터)하는 경우는 조기 종료합니다(수치적 안전장치).

---

## 7) 유효성(Validity)

- 점/접선/보조계수 \((\mathbf{P}_1,\mathbf{D}_1,\mathbf{P}_2,\mathbf{D}_2,\mathbf{C},\mathbf{D})\)가 유효하고, 차원 \(d\in\{2,3\}\)이면 유효.
- 오프셋 생성자에서는 분모/외적 크기 등의 수치적 퇴화를 방지하는 검사 포함.

---

## 8) 구현 포인트 요약

- Hermite ↔ Bezier 변환을 적극 활용: 평가/바운딩/NURBS 변환이 간단해짐
- 미리 저장한 \(\mathbf{C},\mathbf{D}\)로 빠른 평가/도함수 계산
- 평면 오프셋은 제어다각형 기반의 기하식으로 세 가지 경우 분기
- 자연 구간 \([0,1]\), 필요 시 `ON_Xform`으로 변환 후 AABB 갱신

---

## 9) 참고: 간단한 사용 예

```cpp
// 3D Hermite 만들기
ON_3dPoint  P1(0,0,0), P2(1,1,0);
ON_3dVector D1(1,0,0), D2(0,1,0);
ON_HermiteCurve H(P1,D1,P2,D2,3);

// 평가 (점과 1~3차 도함수)
ON_3dPoint d[4];
H.Evaluate(0.25, 4, true, d); // d[0]=pos, d[1]=1st, d[2]=2nd, d[3]=3rd

// Bezier/NURBS 변환
ON_BezierCurve B;
ON_NurbsCurve  N;
H.GetBezierCurve(B);
H.GetNurbForm(N);

// 평면 오프셋 (법선 n, 거리 delta)
ON_3dVector n(0,0,1);
double delta = 0.1;
ON_HermiteCurve Hoff(3, H, n, delta);
```

---

**GitHub 수식 팁**  
- 인라인 수식: `\( a+b \)`  
- 블록 수식: 빈 줄 + `$$`로 감싸기  
- 파이프(`|`)·역슬래시(`\`)가 많은 식은 코드블록 대신 수식 블록을 권장
