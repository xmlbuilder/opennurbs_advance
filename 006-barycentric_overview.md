# Barycentric Coordinates

## 📖 개요
Barycentric Coordinates는 **삼각형 내부의 점을 표현하는 좌표계**입니다.  
삼각형의 정점을 (1,0,0), (0,1,0), (0,0,1)로 두었을 때, 삼각형 내부 점의 위치를 세 정점에 떨어진 비율로 나타냅니다.

- 접두사 `bary-` = weight (가중치)
- 접미사 `-centric` = 중심  
즉, 질량 중심을 의미하는 좌표계이지만 **“질량” 그 자체와는 무관**합니다.

---

## 🧮 수학적 정의
점 P를 삼각형의 세 정점 A, B, C의 선형 조합으로 표현합니다.

$$
P=\alpha A+\beta B+\gamma C,\quad \alpha+\beta+\gamma=1
$$

여기서 α, β, γ는 삼각형의 영역 비율을 나타냅니다.

$$
\alpha = \frac{A_A}{A_A + A_B + A_C}, \quad
\beta = \frac{A_B}{A_A + A_B + A_C}, \quad
\gamma = \frac{A_C}{A_A + A_B + A_C}
$$

- α, β, γ ≥ 0 일 때, 점 P는 삼각형 내부에 존재합니다.
- α, β, γ의 합은 항상 1입니다.

---

## 📊 예시: 좌표계 비교
일반적으로 많이 쓰는 Cartesian Coordinates(데카르트 좌표계)는 **절대적 위치**를 나타냅니다.  
반면, Barycentric Coordinates는 **세 개의 기준점(A,B,C)에 대한 비율**을 나타냅니다.

---

## 📐 시각 자료

### 기본 삼각형과 내부 점
![bary1](/image/Barycentric1.png)




### 영역 비율로 표현된 공식
![bary2](/image/Barycentric2.png)





### 그래픽스 응용 (RGB, 텍스처 매핑)
![bary3](/image/Barycentric3.png)




---

## 🎨 그래픽스 응용
Barycentric Coordinates는 컴퓨터 그래픽스에서 널리 활용됩니다.

- **RGB 색상 보간** : 삼각형 꼭짓점의 색상을 기준으로 내부 점의 색상 계산
- **Texture Mapping** : Mesh의 polygon에 텍스처 좌표를 매핑
- **Intersection & Collision** : 삼각형 내부에 점이 포함되는지 판정

---

## 📝 요약
- Barycentric Coordinates = 삼각형 내부 좌표 표현 방식
- 특징: (α, β, γ)의 합 = 1
- 그래픽스/로보틱스/FEM에서 중요하게 사용됨



---

# 🔧 C++ / OpenGL 예제

아래 예제는 **삼각형 내부에서 Barycentric 보간으로 색상/텍스처를 계산**하는 두 가지 방법을 보여줍니다.

## 1) GLSL에서 에지 함수로 Barycentric 계산 (프래그먼트 셰이더)

> OpenGL의 기본 보간기를 쓰지 않고, **직접** 무게(α,β,γ)를 계산해 색을 만들기.

### Vertex Shader (`bary.vert`)
```glsl
#version 330 core
layout(location = 0) in vec3 inPos;
layout(location = 1) in vec2 inUV;

out vec2 vUV;
void main() {
  vUV = inUV;
  gl_Position = vec4(inPos, 1.0);
}
```

### Fragment Shader (`bary.frag`)
```glsl
#version 330 core
in vec2 vUV;
out vec4 FragColor;

// 화면공간 삼각형 좌표 (NDC가 아닌 경우엔 동일한 공간으로 사전 변환 필요)
uniform vec2 A;  // vertex A xy
uniform vec2 B;  // vertex B xy
uniform vec2 C;  // vertex C xy

// 현재 프래그먼트의 화면 좌표 (픽셀 좌표 또는 NDC에 맞춰 A,B,C와 동일한 좌표계여야 함)
uniform vec2 P;

vec3 barycentric(vec2 a, vec2 b, vec2 c, vec2 p)
{
    // 에지 함수 기반 (signed area)
    float area = (b.x - a.x)*(c.y - a.y) - (b.y - a.y)*(c.x - a.x);
    float a0   = (b.x - p.x)*(c.y - p.y) - (b.y - p.y)*(c.x - p.x);
    float a1   = (c.x - p.x)*(a.y - p.y) - (c.y - p.y)*(a.x - p.x);
    float a2   = (a.x - p.x)*(b.y - p.y) - (a.y - p.y)*(b.x - p.x);
    // α, β, γ (합이 1)
    return vec3(a0, a1, a2) / area;
}

void main()
{
    vec3 w = barycentric(A, B, C, P); // (α, β, γ)
    // 삼각형 내부 판정: 모두 [0,1] 근처
    if (any(lessThan(w, vec3(0.0))) || any(greaterThan(w, vec3(1.0)))) discard;

    // 꼭짓점 색상 (A=red, B=green, C=blue) 가정
    vec3 Ca = vec3(1,0,0);
    vec3 Cb = vec3(0,1,0);
    vec3 Cc = vec3(0,0,1);

    vec3 color = w.x * Ca + w.y * Cb + w.z * Cc;
    FragColor = vec4(color, 1.0);
}
```

> 실사용 시에는 `P`를 `gl_FragCoord.xy`에서 뽑고, `A/B/C`도 같은 좌표계(윈도우 좌표)로 넘겨주면 됩니다.

---

## 2) C++에서 Barycentric으로 텍스처 좌표 보간 (CPU)

> CPU에서 임의의 점 `P`가 삼각형 `(A,B,C)` 안에 있는지 검사하고, `uv`를 보간합니다.

```cpp
#include <array>
#include <cmath>
#include <optional>

struct Vec2 { float x, y; };
struct Vec3 { float x, y, z; };

static float edgeFunction(const Vec2& a, const Vec2& b, const Vec2& p) {
    return (b.x - a.x)*(p.y - a.y) - (b.y - a.y)*(p.x - a.x);
}

struct BaryResult {
    float a, b, c; // α, β, γ
};

std::optional<BaryResult> barycentric(const Vec2& A, const Vec2& B, const Vec2& C, const Vec2& P) {
    float area = edgeFunction(A, B, C);
    if (std::abs(area) < 1e-12f) return std::nullopt; // degenerate

    float a = edgeFunction(B, C, P) / area; // α
    float b = edgeFunction(C, A, P) / area; // β
    float c = 1.0f - a - b;                 // γ
    if (a < 0.f || b < 0.f || c < 0.f) return std::nullopt; // outside
    return BaryResult{a,b,c};
}

// UV 보간
Vec2 interpolateUV(const Vec2& uvA, const Vec2& uvB, const Vec2& uvC, const BaryResult& w) {
    return { w.a*uvA.x + w.b*uvB.x + w.c*uvC.x,
             w.a*uvA.y + w.b*uvB.y + w.c*uvC.y };
}
```

---

## 3) 최소 C++/OpenGL 파이프라인 스니펫 (GLFW + GLAD)

> 삼각형을 그리고, 꼭짓점 색을 자동 보간(하드웨어)하거나 위 셰이더로 직접 보간할 수 있습니다.

```cpp
// CMake: find_package(glfw3 REQUIRED) 등 환경 구성 필요
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <vector>

int main() {
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    GLFWwindow* win = glfwCreateWindow(800, 600, "Barycentric Demo", nullptr, nullptr);
    glfwMakeContextCurrent(win);
    gladLoadGLLoader((GLADloadproc)glfwGetProcAddress);

    float verts[] = {
        //  x     y     z      u   v
        -0.8f, -0.6f, 0.0f,   0.0f, 0.0f, // A
         0.0f,  0.8f, 0.0f,   0.5f, 1.0f, // B
         0.8f, -0.6f, 0.0f,   1.0f, 0.0f  // C
    };
    unsigned vao, vbo;
    glGenVertexArrays(1, &vao);
    glGenBuffers(1, &vbo);
    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(verts), verts, GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 5*sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 5*sizeof(float), (void*)(3*sizeof(float)));
    glEnableVertexAttribArray(1);

    // 셰이더 컴파일/링크는 생략. 위의 bary.vert / bary.frag 사용.
    // uniform A,B,C,P 설정 필요.

    while (!glfwWindowShouldClose(win)) {
        glClear(GL_COLOR_BUFFER_BIT);
        glBindVertexArray(vao);
        glDrawArrays(GL_TRIANGLES, 0, 3);
        glfwSwapBuffers(win);
        glfwPollEvents();
    }
    glfwTerminate();
    return 0;
}
```

> 텍스처 매핑을 하려면 `sampler2D`로 텍스처를 바인딩하고, `interpolateUV` 로 CPU에서 샘플링하거나, 일반적으로는 **GPU에서 자동 보간된 `vUV`로 샘플링**합니다.

---

## ✅ 팁
- 수치 안정성을 위해 면적/에지 함수 계산 시 `double` 사용을 고려하세요.
- 스크린 픽셀 좌표계를 사용할 때는 **A,B,C,P가 모두 동일한 좌표계**인지 확인하세요.
- 뒤집힌 삼각형(시계/반시계)에 따라 부호가 바뀔 수 있으니 절댓값과 부호 처리를 주의하세요.
