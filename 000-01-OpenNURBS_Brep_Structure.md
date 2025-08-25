# 제 5 장. OpenNURBS의 B-rep 구조

## 5.1 B-rep 개요
Boundary Representation (B-rep)은 3차원 형상을 **경계 곡선과 곡면**으로 표현하는 방식이다.  
CAD 커널에서는 B-rep을 통해 기하학적 형상(Geometry)과 위상학적 구조(Topology)를 명확히 구분한다.

- **Geometry**: 실제 공간의 곡선과 곡면 (NURBS 기반 데이터)
- **Topology**: Geometry 간의 연결과 경계를 정의하는 구조

---

## 5.2 기하학 (Geometry)

B-rep이 참조하는 기하학 객체들은 다음과 같다:

- `ON_CurveArray m_C2` : 2D Parameter 공간에서의 Trimming 곡선  
- `ON_CurveArray m_C3` : 3D 공간의 Edge 곡선  
- `ON_SurfaceArray m_S` : Face를 정의하는 NURBS 곡면  

즉, Surface 위의 Loop는 2D trimming curve (`m_C2`)로 정의되고,  
Edge는 3D 곡선 (`m_C3`)로 정의된다.

---

## 5.3 위상 (Topology)

Topology는 Geometry의 관계를 조직화한 것이다.  
OpenNURBS에서는 아래의 다섯 가지 요소로 위상이 정의된다:

- **Vertex (꼭짓점)** : Edge의 끝점 (3D 좌표)
- **Edge (간선)** : Vertex를 연결하는 3D 곡선 (`m_C3` 참조)
- **Trim (트림)** : Edge의 파라미터 공간 표현 (2D trimming curve, `m_C2` 참조)
- **Loop (루프)** : 여러 Trim이 모여 형성하는 닫힌 경계 (외곽 루프, 홀 루프)
- **Face (면)** : Loop로 경계가 정의된 Surface (`m_S` 참조)

---

## 5.4 위상 구조 다이어그램

![Brep Structure](/image/brep_structure.png)

- **Face**는 하나 이상의 **Loop**를 가진다.  
- **Loop**는 여러 **Trim**으로 구성된다.  
- **Trim**은 **Edge**와 연결되고, Edge는 **Vertex**를 참조한다.  
- 각 Face는 특정 **Surface**를 참조한다.  

즉, Face ↔ Loop ↔ Trim ↔ Edge ↔ Vertex 구조가 **양방향 참조**로 이어진다.

---

## 5.5 수학적 의미

B-rep 구조의 핵심은 **Surface의 부분 집합을 정의**하는 것이다.  

수학적으로 Face는 다음과 같이 표현할 수 있다:

$$
F = \{ (u,v) \in D \subset \mathbb{R}^2 \mid \gamma_i(u,v) \geq 0, \; i=1,...,n \}
$$

여기서,
- \( D \) : Surface의 전체 파라미터 영역  
- $\( T_i \)$ : Trim 곡선 (Loop를 이루는 곡선들)  
- \( F \) : Trimmed Surface  

즉, Face는 NURBS Surface 위에서 Trim 곡선으로 경계가 잘린 부분이다.

---

## 5.6 예시: Trimmed Surface

### 1단계. Surface 생성
```cpp
ON_NurbsSurface srf = ...; // 기본 NURBS 곡면
```

### 2단계. Loop 정의
```cpp
ON_BrepLoop* loop = brep.NewLoop(ON_BrepLoop::outer, face);
```

### 3단계. Trim 곡선 추가
```cpp
ON_BrepTrim* trim = brep.NewTrim(edge, loop, dir);
```

### 4단계. Edge / Vertex 생성
```cpp
ON_BrepVertex& v0 = brep.NewVertex(pt0, tol);
ON_BrepVertex& v1 = brep.NewVertex(pt1, tol);
```

### 5단계. Face 확정
```cpp
ON_BrepFace& face = brep.NewFace(srf);
```

---

## 5.7 요약
- **Geometry**는 NURBS 곡선/곡면 데이터  
- **Topology**는 Geometry 간 연결 구조 (Face, Loop, Trim, Edge, Vertex)  
- Trimmed Surface는 **Surface + Loop(Trim curves)** 로 정의된다.  
