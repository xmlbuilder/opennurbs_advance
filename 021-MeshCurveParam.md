# `mesh_curve_param` 알고리즘 기술 문서

## 개요
`mesh_curve_param` 은 1차원 곡선(예: B-Rep Edge, NURBS Curve)을 **선분 세그먼트(polyline)**로 적응적으로 분할하는 함수이다.  
이 알고리즘은 곡선 기반 메시 생성(1D → 2D → 3D 확장)에서 가장 기초적인 단계로, 주어진 곡선의 파라미터 구간 `[u0, u1]`을 따라 **Chordal Error**, **Target Length**, **Gradation** 조건을 만족하는 **세그먼트 분할**을 수행한다.

---

## 입력/출력 시그니처

```cpp
template <class Curve>
int mesh_curve_param(
   const Curve& C,
   DoubleMat3& pos3D, unsigned N0, unsigned N1, double u0, double u1,
   double& h0, double& h1, double target_h,
   bool force_even, unsigned min_n, unsigned max_n,
   double max_chordal_error, double min_h, unsigned chordal_control_type,
   double max_gradation,
   UIntMatE& connectE, DoubleVec& pos1D,
   unsigned high_order_type = 0, unsigned max_bgm_remeshings = 4);
```

## 주요 파라미터
- Curve& C : 곡선을 평가할 수 있는 객체 (get_3D_coordinates, get_local_bases, get_local_curvatures 제공).
- pos3D : 3D 점들을 누적 저장하는 컨테이너. (출력 포함)
- N0, N1 : 시작/끝점 인덱스(pos3D 내).
- u0, u1 : 시작/끝 파라미터 (정규화 혹은 실제 도메인).
- h0, h1 : 시작/끝 세그먼트 길이 (출력: 실제 세그 길이 기반 업데이트).
- target_h : 목표 세그먼트 길이.
- min_h : 최소 허용 세그 길이.
- max_chordal_error : 최대 허용 현오차.
- chordal_control_type : 분할 기준 선택
    - 0 : 길이만
    - 1 : 현오차만
    - 2 : 길이 + 현오차 (AND)
- max_gradation : 인접 세그 길이비 허용 한계 (e.g. 1.6).
- connectE : 세그먼트 연결 결과 (E2).
- pos1D : 샘플링된 파라미터 값들 저장.
- high_order_type : 고차 재투영 옵션 (예: 2차/3차 보간).
- max_bgm_remeshings : BGM 유사 반복 횟수.

## 알고리즘 흐름

```mermaid
flowchart TD
  A[입력: 곡선 C, u0,u1, target_h, max_ce, gradation] --> B[초기 세그 생성 (끝점 2개)]
  B --> C[우선순위 큐 초기화 (세그 score)]
  C --> D{큐 empty?}
  D -- 아니오 --> E[큐 top 세그 꺼냄]
  E --> F{길이/현오차 OK?}
  F -- 예 --> G[큐 유지 (다음 세그 검사)]
  F -- 아니오 --> H[중점 분할 → pos3D, pos1D 추가]
  H --> I[좌/우 세그 score 계산 후 큐 삽입]
  I --> D
  D -- 예 --> J[세그 집합 정렬]
  J --> K[그라데이션 검사 + 필요 시 추가 분할]
  K --> L[최종 연결 connectE 생성, h0/h1 업데이트]
  L --> M[출력 완료]
```

## 핵심 아이디어

### 1. Chordal Error 기반 분할
    - 현오차 = 곡선 중점과 직선 중점의 차이
    - chordal_error > max_chordal_error → 분할 필요

### 2. Target Length 기반 분할
    - 세그 길이 > target_h → 분할 필요
    - min_h 보다 짧아지면 분할 중단

### 3. Best-First Refinement (우선순위 큐)
    - 세그 점수 score = max(L/target_h, chordal/max_ce)
    - 점수가 가장 큰 세그먼트부터 분할 → 빠른 수렴
### 4. Gradation 제한 (BGM 스타일)
    - 인접 세그 길이비 > max_gradation → 중점 삽입
    - 2~3회 스윕으로 전체 길이 분포 균질화
### 5. High-Order 보간(옵션)
    - high_order_type>0 이면 중간노드를 다시 곡선에 투영
    - 더 정확한 분할 가능 (현재는 기본 noop)


## 품질 지표
- Chordal Error : 중점 오차 (직선 근사 정확도)
- 세그 길이 분포 : 평균 길이, 최대 길이, min/max ratio
- Gradation : 인접 세그 길이비 (예: ≤1.6 권장)
- 노드 수 : target_h 대비 노드 수 적정성


## 사용 예시
```cpp
// OpenNURBS Edge Curve 예시
const ON_Curve* occ = brep.EdgeCurve(ei);
OpenNurbsCurveGen C(occ, /*normalized=*/true);

DoubleMat3 pos3D;
UIntMatE connectE;
DoubleVec pos1D;

// 시작/끝점 추가
unsigned N0 = (unsigned)pos3D.size(); pos3D.push_back({A.x,A.y,A.z});
unsigned N1 = (unsigned)pos3D.size(); pos3D.push_back({B.x,B.y,B.z});

double u0=0.0, u1=1.0;
double h0=0, h1=0;

mesh_curve_param(C, pos3D, N0, N1, u0, u1,
                 h0, h1, 0.2,
                 false, 2, 100000,
                 0.01, 0.02, 2,
                 1.6, connectE, pos1D,
                 0, 4);

```
## 출력:

- pos3D: 분할된 3D 노드 좌표
- connectE: 인덱스 연결 정보
- pos1D: 샘플 파라미터 값들

## 요약
- mesh_curve_param 은 Chordal Error + 길이 제약 + Gradation 을 동시에 만족하는 적응적 1D 메싱 루틴이다.
- 우선순위 큐 기반 분할로 빠르게 수렴하며, BGM 유사 반복으로 세그 길이 분포를 보정한다.
- OpenNURBS와 같은 CAD 커널의 ON_Curve 를 어댑터(OpenNurbsCurveGen)로 연결하면 바로 곡선 메싱이 가능하다.

## 개념도



## 전체 코드
```cpp
using DoubleVec = std::vector<double>;
using DoubleMat3 = std::vector<std::array<double, 3>>;
using DoubleMat2 = std::vector<std::array<double, 2>>;
using UIntMatE = std::vector<std::array<unsigned, 2>>;
using TriMat = std::vector<std::array<unsigned, 3>>;


static inline double len3(const std::array<double, 3>& p) {
  return std::sqrt((std::max)(0.0, p[0] * p[0] + p[1] * p[1] + p[2] * p[2]));
}

static inline std::array<double, 3> sub3(const std::array<double, 3>& a, const std::array<double, 3>& b)
{
  return { a[0] - b[0],a[1] - b[1],a[2] - b[2] };
}

static inline std::array<double, 3> add3(const std::array<double, 3>& a, const std::array<double, 3>& b)
{
  return { a[0] + b[0],a[1] + b[1],a[2] + b[2] };
}

static inline std::array<double, 3> mul3(const std::array<double, 3>& a, double s)
{
  return { a[0] * s,a[1] * s,a[2] * s };
}
static inline double dot3(const std::array<double, 3>& a, const std::array<double, 3>& b)
{
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static inline std::array<double, 3> cross3(const std::array<double, 3>& a, const std::array<double, 3>& b)
{
  return { a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0] };
}

static inline std::array<double, 3> ONP(const ON_3dPoint& P) 
{ 
  return { P.x, P.y, P.z }; 
}

static inline double norm3(const ON_3dVector& v) 
{ 
  return std::sqrt((std::max)(0.0, v.x * v.x + v.y * v.y + v.z * v.z)); 
}
struct Seg {
  double u0, u1;       // 파라미터 구간
  unsigned i0, i1;     // pos3D 인덱스
  double chord_err;    // 중점 chordal error
  double length;       // |P1-P0|
  double score;        // refine 우선순위 (높을수록 먼저 split)
  int    gen;          // 세대(디버그/추적용)
};

template <class Curve>
static inline void eval_point(const Curve& C, double u, std::array<double, 3>& P)
{
  DoubleVec U{ u };
  DoubleMat3 X; X.reserve(1);
  C.get_3D_coordinates(U, X);
  P = X[0];
}

template <class Curve>
static inline double compute_chordal_mid(const Curve& C, double u0, double u1,
  const std::array<double, 3>& P0,
  const std::array<double, 3>& P1)
{
  const double um = 0.5 * (u0 + u1);
  std::array<double, 3> Pm; eval_point(C, um, Pm);
  auto Lm = mul3(add3(P0, P1), 0.5);
  return len3(sub3(Pm, Lm));
}

static inline double safe_div(double a, double b) {
  return (b > 0.0) ? (a / b) : (a > 0.0 ? std::numeric_limits<double>::infinity() : 0.0);
}

template <class Curve>
int mesh_curve_param(
  const Curve& C,
  std::vector<std::array<double, 3>>& pos3D,
  unsigned N0,
  unsigned N1,
  double u0,
  double u1,
  double& h0,
  double& h1,
  double target_h,
  bool force_even,
  unsigned min_n,
  unsigned max_n,
  double max_chordal_error,
  double min_h,
  unsigned chordal_control_type,
  double max_gradation,
  std::vector<std::array<unsigned, 2>>& connectE,
  std::vector<double>& pos1D,
  unsigned high_order_type = 0,
  unsigned max_bgm_remeshings = 4)
{
  assert(N0 < pos3D.size() && N1 < pos3D.size());
  if (N0 >= pos3D.size() || N1 >= pos3D.size()) return -1;

  // 초기 세그먼트 1개로 시작
  auto P0 = pos3D[N0];
  auto P1 = pos3D[N1];
  double L01 = len3(sub3(P1, P0));
  double ce01 = compute_chordal_mid(C, u0, u1, P0, P1);

  auto score_of = [&](double L, double CE)->double {
    // chordal_control_type:
    // 0: length only, 1: chordal only, 2(default): both의 max
    if (chordal_control_type == 0) return safe_div(L, target_h);
    if (chordal_control_type == 1) return safe_div(CE, max_chordal_error);
    return (std::max)(safe_div(L, (std::max)(min_h, target_h)),
      safe_div(CE, (std::max)(1e-16, max_chordal_error)));
    };

  auto make_seg = [&](double a, double b, unsigned i, unsigned j, int gen)->Seg {
    const auto& A = pos3D[i];
    const auto& B = pos3D[j];
    const double L = len3(sub3(B, A));
    const double CE = compute_chordal_mid(C, a, b, A, B);
    Seg s{ a,b,i,j, CE, L, score_of(L,CE), gen };
    return s;
    };

  struct Cmp { bool operator()(const Seg& x, const Seg& y) const { return x.score < y.score; } };
  std::priority_queue<Seg, std::vector<Seg>, Cmp> pq;
  pq.push(make_seg(u0, u1, N0, N1, 0));

  // 스플릿 루프 (best-first)
  unsigned max_total_nodes = std::max<unsigned>(N1 + 1, (unsigned)pos3D.size());
  max_total_nodes = std::min<unsigned>(max_total_nodes + max_n, 5u * max_n + 1024u);

  unsigned splits = 0;
  const unsigned max_splits = 100000u; // 안전장치
  while (!pq.empty())
  {
    Seg s = pq.top();
    // 충분히 만족 → 종료
    bool chord_ok = (max_chordal_error <= 0.0) ? true : (s.chord_err <= max_chordal_error);
    bool len_ok = (s.length <= (std::max)(min_h, target_h));
    bool ok = false;
    if (chordal_control_type == 0) ok = len_ok;
    else if (chordal_control_type == 1) ok = chord_ok;
    else ok = (chord_ok && len_ok);

    if (ok) break; // 최상위도 만족 → 전체 종료

    pq.pop();

    // min_h 이하로 너무 짧아지면 더 이상 split 금지
    if (s.length <= (std::max)(1e-12, min_h)) continue;

    // 중점 추가
    double um = 0.5 * (s.u0 + s.u1);
    std::array<double, 3> Pm; eval_point(C, um, Pm);
    pos3D.push_back(Pm);
    unsigned im = (unsigned)pos3D.size() - 1;
    pos1D.push_back(um);

    // 좌/우 세그 다시 계산해서 push
    pq.push(make_seg(s.u0, um, s.i0, im, s.gen + 1));
    pq.push(make_seg(um, s.u1, im, s.i1, s.gen + 1));

    if (++splits > max_splits) break;
    if (pos3D.size() >= max_total_nodes) break;
  }

  // 우선순위 큐의 모든 세그를 꺼내서 파라미터 순으로 정렬 → 최종 polyline
  std::vector<Seg> segs;
  segs.reserve(pq.size());
  while (!pq.empty()) { segs.push_back(pq.top()); pq.pop(); }
  if (segs.empty()) {
    // 최소 보장
    segs.push_back(make_seg(u0, u1, N0, N1, 0));
  }
  std::sort(segs.begin(), segs.end(), [](const Seg& a, const Seg& b) { return a.u0 < b.u0; });

  // 인접 세그 그라데이션 제한 (BGM 스타일로 2~3회 스윕)
  auto split_mid = [&](const Seg& s)->std::pair<Seg, Seg> {
    double um = 0.5 * (s.u0 + s.u1);
    std::array<double, 3> Pm; eval_point(C, um, Pm);
    pos3D.push_back(Pm);
    unsigned im = (unsigned)pos3D.size() - 1;
    pos1D.push_back(um);
    return { make_seg(s.u0,um,s.i0,im,s.gen + 1), make_seg(um,s.u1,im,s.i1,s.gen + 1) };
    };

  unsigned grad_passes = std::min<unsigned>(max_bgm_remeshings, 4u);
  for (unsigned it = 0; it < grad_passes; ++it) {
    bool changed = false;
    std::vector<Seg> next; next.reserve(segs.size() * 2);
    for (size_t k = 0; k < segs.size(); ++k) {
      const Seg& s = segs[k];
      double Ls = s.length;
      double Lprev = (k > 0) ? segs[k - 1].length : Ls;
      double Lnext = (k + 1 < segs.size()) ? segs[k + 1].length : Ls;

      double r1 = (Lprev > 0) ? (std::max)(Ls / Lprev, Lprev / Ls) : 1.0;
      double r2 = (Lnext > 0) ? (std::max)(Ls / Lnext, Lnext / Ls) : 1.0;
      double worst = (std::max)(r1, r2);

      if (worst > (std::max)(1.0, max_gradation) && Ls > (std::max)(1e-12, min_h)) {
        auto p = split_mid(s);
        next.push_back(p.first);
        next.push_back(p.second);
        changed = true;
      }
      else {
        next.push_back(s);
      }
    }
    if (!changed) break;
    // 다음 패스 위해 재정렬(파라미터 순)
    std::sort(next.begin(), next.end(), [](const Seg& a, const Seg& b) { return a.u0 < b.u0; });
    segs.swap(next);
  }

  // 출력: connectE(연결) & h0/h1
  connectE.clear();
  if (segs.size() == 1) {
    connectE.push_back({ segs[0].i0, segs[0].i1 });
  }
  else {
    // 세그들을 이어붙인 polyline으로 변환 (오버랩/틈 없이 정렬되어 있음)
    connectE.push_back({ segs.front().i0, segs.front().i1 });
    for (size_t k = 1; k < segs.size(); ++k) {
      // 앞 세그의 i1 == 뒤 세그의 i0 인 상태가 대부분이지만,
      // 안전하게 연결관계를 명시
      if (segs[k - 1].i1 != segs[k].i0) {
        // 중간이 어긋났다면 보정 연결 (이 케이스는 거의 없음)
        connectE.push_back({ segs[k - 1].i1, segs[k].i0 });
      }
      connectE.push_back({ segs[k].i0, segs[k].i1 });
    }
  }

  // h0/h1(끝단 길이 메트릭) 업데이트
  if (!segs.empty()) {
    double L0 = segs.front().length;
    double L1 = segs.back().length;
    h0 = (std::max)(min_h, (std::min)(L0, target_h));
    h1 = (std::max)(min_h, (std::min)(L1, target_h));
  }
  return 0;
}


```
