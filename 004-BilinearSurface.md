# 🌀 ON_MakeBilinearSurface Example

이 문서는 OpenNURBS를 이용해 **Bilinear NURBS Surface**를 생성하고, 이를
메쉬로 변환 후 BDF 형식으로 내보내는 예제를 설명합니다.

## 📐 Bilinear Surface 생성 함수

``` cpp
bool ON_MakeBilinearSurface(
  const ON_3dPoint& P00, 
  const ON_3dPoint& P10,
  const ON_3dPoint& P01, 
  const ON_3dPoint& P11,
  ON_NurbsSurface& srf
)
{
  if (!srf.Create(3, false, 1+1, 1+1, 2, 2)) return false;

  // u, v 모두 차수 1, clamped knots
  srf.SetKnot(0, 0, 0.0); srf.SetKnot(0, 1, 1.0);
  srf.SetKnot(1, 0, 0.0); srf.SetKnot(1, 1, 1.0);

  // CV 배치 (i=u index, j=v index)
  srf.SetCV(0, 0, P00); srf.SetCV(1, 0, P10);
  srf.SetCV(0, 1, P01); srf.SetCV(1, 1, P11);

  return srf.IsValid();
}
```

-   차수: U, V 모두 1 (Linear)
-   노트: \[0, 1\] 구간 Clamped Knot
-   Control Points:
    -   (0,0) = P00, (1,0) = P10
    -   (0,1) = P01, (1,1) = P11

------------------------------------------------------------------------

## 🚀 실행 예제

``` cpp
int main(int argc, const char* argv[])
{
  ON::Begin();
  ON_TextLog log;

  ON_NurbsSurface srf;
  ON_MakeBilinearSurface({ 0,0,0 }, { 10,0,2 }, { 0,8,1 }, { 10,8,0 }, srf);

  srf.Dump(log);
  srf.IsValid(&log);

  ON_Mesh mesh;
  ON_MeshParameters mp;
  mp.SetGridMaxCount(20);
  srf.CreateMesh(mp, &mesh);
  ON_ExportBDF bdf;
  bdf.Run(L"D:\\Temp\\Test1.bdf", &mesh);
 
  ON::End();
  return 0;
}
```

------------------------------------------------------------------------

## ✅ 실행 결과

    ON_NurbsSurface dim = 3 is_rat = 0
            order = 2 X 2 cv_count = 2 X 2
    Knot Vector 0 ( 2 knots )
    index                     value  mult       delta
        0                        0     1
        1                        1     1           1
    Knot Vector 1 ( 2 knots )
    index                     value  mult       delta
        0                        0     1
        1                        1     1           1
    Control Points  4 non-rational points
      index               value
      CV[ 0][ 0] (0, 0, 0)
      CV[ 0][ 1] (0, 8, 1)
      CV[ 1][ 0] (10, 0, 2)
      CV[ 1][ 1] (10, 8, 0)

## 참고 그림
![bilinear surface](/image/bilinear_surface.jpg)

------------------------------------------------------------------------

## 🔧 유틸리티 생성자

``` cpp
ON_3dPoint(std::initializer_list<double> vec);
ON_3dPoint::ON_3dPoint(std::initializer_list<double> vec)
{
  if (vec.size() == 3)
  {
    int i = 0;
    for (auto itr = vec.begin(); itr != vec.end(); itr++)
    {
      if (i == 0) x = (*itr);
      else if (i == 1) y = (*itr);
      else if (i == 2) z = (*itr);
      i++;
    }
  }
}
```

-   `std::initializer_list`를 활용하여 `{x, y, z}` 형태로 손쉽게
    `ON_3dPoint`를 생성 가능.

------------------------------------------------------------------------

## 📂 요약

-   `ON_MakeBilinearSurface()` : 4개의 점으로 Bilinear Surface 생성
-   `CreateMesh()` : Surface → Mesh 변환 (OpenNurbs에 없는 함수 추후에 공개)
-   `ON_ExportBDF` : 메쉬를 BDF 형식으로 내보내기(OpenNurbs에 없는 함수 추후에 공개)
-   `ON_3dPoint({x,y,z})` : 초기화 리스트 생성자 제공

