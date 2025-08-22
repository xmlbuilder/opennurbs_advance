# 🧩 OpenNURBS Advance

**목표:** OpenNURBS API를 기반으로 Rhino가 제공하는 **고급 NURBS 모델링
기능**(스윙/로프트/스윕/쌍선형/병진 스윕 등)을\
"The NURBS Book"의 이론에 따라 **재구현/확장**하는 실험 프로젝트입니다.\
학습/연구 목적의 예제와 구현 노트, 그리고 작은 데모 코드를 함께
제공합니다.\
(레포 소개: "OpenNURBS API 활용, The NURBS Book 이론 기반" 문구 참조)

------------------------------------------------------------------------

## 📂 구성

문서와 코드는 주제별 파일로 정리됩니다. (일부 예)

-   `001.swung_surface.md` --- **Swung Surface** 이론/알고리즘/예제
-   `002.loft_surface.md` --- **Loft Surface** 생성 전략과 파라미터 선택
-   `003-01.sweep_surface.md` --- **Sweep Surface** (1/2):
    개념/조건/경계
-   `003-02-sweep_surface_demo.cpp` --- **Sweep 데모 코드** (C++)
-   `004-BilinearSurface.md` --- **Bilinear Surface** (쌍선형 보간)
-   `005.MakeTranslationalSweep.md` --- **Translational Sweep** (병진
    스윕)

> 파일 목록은 레포 루트에서 확인할 수 있습니다.

------------------------------------------------------------------------

## 🛠️ 빌드 & 실행 (예시)

> **전제:** OpenNURBS SDK가 준비되어 있고(`opennurbs.h` 등), C++17 이상
> 컴파일러를 사용합니다.\
> OS/툴체인은 자유롭게 선택하세요(Visual Studio / Clang / GCC / CMake
> 등).

### 1) 단일 파일 데모(예: 스윕 데모)

``` bash
# 예) GCC/Clang (경로/라이브러리는 환경에 맞게 수정)
g++ -std=c++17 003-02-sweep_surface_demo.cpp -I/path/to/opennurbs/include     -L/path/to/opennurbs/lib -lopennurbs -o sweep_demo

./sweep_demo
```

### 2) CMake 사용(권장)

``` cmake
# CMakeLists.txt (예시)
cmake_minimum_required(VERSION 3.20)
project(opennurbs_advance CXX)

set(CMAKE_CXX_STANDARD 17)
include_directories(/path/to/opennurbs/include)
link_directories(/path/to/opennurbs/lib)

add_executable(sweep_demo 003-02-sweep_surface_demo.cpp)
target_link_libraries(sweep_demo opennurbs)
```

``` bash
mkdir build && cd build
cmake ..
cmake --build . --config Release
./sweep_demo
```

------------------------------------------------------------------------

## 📖 학습 가이드

각 문서는 다음 흐름으로 작성됩니다.

1.  **이론 개요**: 정의, 파라미터화, 경계 조건, 제어점/노트 벡터 등 핵심
    수식 정리\
2.  **알고리즘/절차**: 입력 곡선/단면/레일, 샘플링/병합,
    표준화(normalization)\
3.  **코드 스니펫**: OpenNURBS API를 사용하는 최소 예제\
4.  **주의/함정**: 파라미터 정렬, 곡률/연속성(C0/C1/C2), 수치 안정성,
    예외 케이스\
5.  **확장 아이디어**: Rhino 기능과의 차이, 성능/품질 옵션,
    후처리(삼각분할/메쉬화) 등

------------------------------------------------------------------------

## 🔎 로드맵

-   [x] Swung / Loft / Sweep / Bilinear / Translational Sweep 문서 초안
-   [ ] Sweep (2/2): **Rail/Section** 변형·Twist·Scaling 옵션
-   [ ] Loft 고급: Matching/Refit, 클로즈드 섹션 처리, 연속성
    제어(G0/C1/C2)
-   [ ] Swung 고급: 회전체 일반화, Param coupling 전략
-   [ ] Nurbs 기본기 보강: Knot/Multiplicity, De Boor, Global/Local
    Interp
-   [ ] 메쉬화: Adaptive tessellation, curvature-based sampling
-   [ ] 품질/성능: 톨러런스/유효 자릿수, Degenerate edge/trim 보정
-   [ ] 테스트 데이터셋 & 시각화 스크립트(이미지/뷰어)

------------------------------------------------------------------------

## 🤝 Contributing

-   이론 보강, 예제 코드(PR), 반례/테스트 데이터 환영합니다.\
-   Issue에 아이디어·질문을 남겨주세요.\
-   스타일: C++17, `clang-format` 권장 (Google/LLVM 스타일 중 택1)

------------------------------------------------------------------------

## 📚 참고

-   **OpenNURBS** API 문서/헤더\
-   **The NURBS Book** (Piegl & Tiller) --- 본 레포의 이론적 기반\
-   Rhino 사용자 문서(기능 레벨 비교 참고용)

------------------------------------------------------------------------

## 📜 License

MIT (문서/코드에 별도 표기가 있는 경우 해당 라이선스 우선)
