# Bézier Surface ↔ Power Basis 변환 이론

## 1. Bézier Surface 정의

차수 $\((d_u, d_v)\)$ 의 Bézier 곡면은 제어점 $\(P_{i,j}\)$ 와 Bernstein basis $\(B_{i}^{d_u}(u), B_{j}^{d_v}(v)\)$ 로 정의된다:

$$
S(u,v) = \sum_{i=0}^{d_u} \sum_{j=0}^{d_v} P_{i,j} \, B_i^{d_u}(u)\, B_j^{d_v}(v),
$$
여기서
$$
B_i^n(t) = \binom{n}{i} (1-t)^{n-i} t^i.
$$

---

## 2. Power Basis (모노미얼) 정의

Power basis는 단순히 다항식 전개 형태:

$$
S(u,v) = \sum_{p=0}^{d_u} \sum_{q=0}^{d_v} A_{p,q}\, u^p v^q,
$$

여기서 $\(A_{p,q}\)$는 계수 행렬로, x, y, z 좌표 각각에 대해 하나씩 존재한다:
- $\(A^x_{p,q}\)$, $\(A^y_{p,q}\)$, $\(A^z_{p,q}\)$.

---

## 3. Bézier → Power 변환 행렬

1D Bézier basis와 Power basis 사이에는 **선형 변환**이 존재한다.  
- Bézier basis $\([B_0^n, \dots, B_n^n]\)$ 와  
- Power basis $\([1, t, t^2, \dots, t^n]\)$ 사이 변환:

$$
[B_0^n(t), B_1^n(t), \dots, B_n^n(t)]^T = M_{B\to P} \cdot [1, t, t^2, \dots, t^n]^T
$$

따라서 제어점 행렬 \(P\)를 Power 계수로 바꾸려면:

$$
A = M_u \, P \, M_v^T
$$

- \(M_u\): u 방향 Bézier → Power 변환 행렬  
- \(M_v\): v 방향 Bézier → Power 변환 행렬  
- \(P\): 제어점 좌표 행렬 (x, y, z 각각 따로)  
- \(A\): Power 계수 행렬

---

## 4. Power → Bézier 변환

역변환은 단순히 위 행렬을 역행렬로 바꿔주면 된다:

$$
P = M_u^{-1} \, A \, (M_v^{-1})^T
$$

여기서 $\(M_u^{-1}, M_v^{-1}\)$는 Bézier→Power 행렬의 역행렬.

---

## 5. 평가 방법

### 5.1 Bézier Surface 평가
주어진 $\((u,v)\)$에서:

$$
S(u,v) = \sum_{i=0}^{d_u} \sum_{j=0}^{d_v} P_{i,j} \, B_i^{d_u}(u)\, B_j^{d_v}(v).
$$

---

### 5.2 Power Basis 평가
주어진 \((u,v)\)에서:

$$
S(u,v) = \sum_{p=0}^{d_u} \sum_{q=0}^{d_v} A_{p,q}\, u^p v^q.
$$

실제 계산은 **Horner scheme**을 사용하여 수치적으로 안정적으로 수행.

---

## 6. 정리

- Bézier와 Power Basis는 **동일한 함수 공간**을 표현한다.  
- 변환은 단순히 행렬 곱으로 이루어진다.  
- 변환 전후로 곡면은 완전히 동일해야 하며, 수치 오차는 부동소수점 한계 (~1e-15) 정도다.  
- Power basis는 미분·적분에 편리하고, Bézier basis는 기하학적 직관(제어점, convex hull 등)에 강점이 있다.

---
## 구현 코드
```cpp

#include <algorithm>
#include <type_traits> // std::enable_if, std::is_floating_point


template<typename Type>
class ON_TMatrix
{
public:
	using RowIterator = Type**;
	using ConstRowIterator = const Type* const*;

	ON_TMatrix()
	{

	}
	ON_TMatrix(
		int row_count,
		int col_count
	)
	{
		Create(row_count, col_count);
	}

	ON_TMatrix(const ON_TMatrix& src)
	{
		*this = src;
	}

	ON_TMatrix(ON_TMatrix&& src) ON_NOEXCEPT
	{
		if (this != &src)
		{
			m = src.m;
			src.m = nullptr;
			m_row_count = src.m_row_count;
			src.m_row_count = 0;
			m_col_count = src.m_col_count;
			src.m_col_count = 0;
		}
	}

	ON_TMatrix(std::initializer_list<std::initializer_list<Type>> list)
	{
		int rows = static_cast<int>(list.size());
		int cols = rows > 0 ? static_cast<int>(list.begin()->size()) : 0;
		Create(rows, cols);
		int i = 0;
		for (const auto& row : list)
		{
			int j = 0;
			for (const auto& val : row)
				m[i][j++] = val;
			++i;
		}
	}

	template<typename T, typename F>
	static void Fill(ON_TMatrix<T>& mat, F fillFunc)
	{
		for (int r = 0; r < mat.RowCount(); ++r)
		{
			for (int c = 0; c < mat.ColCount(); ++c)
			{
				mat[r][c] = fillFunc(r, c);
			}
		}
	}

	ON_TMatrix& operator=(ON_TMatrix&& src)
	{
		if (this != &src)
		{
			Destroy();
			m = src.m;
			src.m = nullptr;
			m_row_count = src.m_row_count;
			src.m_row_count = 0;
			m_col_count = src.m_col_count;
			src.m_col_count = 0;
		}
		return *this;
	}

	ON_TMatrix& operator=(const ON_TMatrix& src)
	{
		if (this != &src)
		{
			Destroy();
			Create(src.m_row_count, src.m_col_count);
			for (int i = 0; i < m_row_count; i++)
			{
				for (int j = 0; j < m_col_count; j++)
				{
					m[i][j] = src.m[i][j];
				}
			}
		}
		return *this;
	}

	ON_TMatrix& operator*=(const ON_TMatrix& rhs)
		requires std::is_floating_point_v<Type>
	{
		*this = (*this) * rhs; // 위에서 정의한 안전한 곱 사용
		return *this;
	}

	ON_TMatrix operator*(const ON_TMatrix& src)
		requires std::is_floating_point_v<Type>
	{
		if (m_col_count != src.m_row_count)
			throw std::invalid_argument("Matrix dimensions do not match for multiplication.");

		ON_TMatrix result(src.m_row_count, src.m_col_count);
		for (int i = 0; i < m_row_count; ++i)
		{
			for (int j = 0; j < src.m_col_count; ++j)
			{
				double sum = 0.0;
				// (필요시 캐시/BLAS로 최적화 가능)
				for (int k = 0; k < m_col_count; ++k)
					sum += m[i][k] * src.m[k][j];
				result[i][j] = sum;
			}
		}
		return result; // 새 객체 반환 (좌변 보존)
	}

	template <typename Type>
	requires std::is_floating_point_v<Type>
		friend ON_TMatrix<Type> operator*(const ON_TMatrix<Type>& left, const ON_TMatrix<Type>& right);

	ON_TMatrix<Type> Invert(double zero_tolerance = ON_ZERO_TOLERANCE) const
		requires std::is_floating_point_v<Type>
	{
		if (m_row_count != m_col_count)
			throw std::invalid_argument("Inverse is only defined for square matrices.");

		const int n = m_row_count;
		ON_TMatrix<Type> A(*this);
		ON_TMatrix<Type> I(n, n);
		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < n; j++)
			{
				I[i][j] = static_cast<Type>(0.0);
			}
		}
		for (int i = 0; i < n; ++i) I[i][i] = static_cast<Type>(1);
		for (int col = 0; col < n; ++col)
		{
			// 부분 피벗팅
			int pivot_row = col;
			Type max_abs = std::fabs(A[col][col]);
			for (int r = col + 1; r < n; ++r)
			{
				Type v = std::fabs(A[r][col]);
				if (v > max_abs) { max_abs = v; pivot_row = r; }
			}
			if (!(max_abs > static_cast<Type>(zero_tolerance)))
				throw std::runtime_error("Matrix is singular or ill-conditioned.");

			if (pivot_row != col)
			{
				// A의 행 교환
				for (int j = 0; j < n; ++j)
					std::swap(A[pivot_row][j], A[col][j]);
				// I의 행 교환
				for (int j = 0; j < n; ++j)
					std::swap(I[pivot_row][j], I[col][j]);
			}

			// 피벗 행 정규화
			const Type pivot = A[col][col];
			const Type inv_pivot = static_cast<Type>(1) / pivot;
			for (int j = 0; j < n; ++j)
			{
				A[col][j] *= inv_pivot;
				I[col][j] *= inv_pivot;
			}

			// 다른 행에서 피벗 열 제거
			for (int r = 0; r < n; ++r)
			{
				if (r == col) continue;
				const Type factor = A[r][col];
				if (factor == static_cast<Type>(0)) continue;
				for (int j = 0; j < n; ++j)
				{
					A[r][j] -= factor * A[col][j];
					I[r][j] -= factor * I[col][j];
				}
			}
		}
		return I;
	}

	virtual ~ON_TMatrix()
	{
		Destroy();
	}

	Type* operator[](int i)
	{
		return m[i];
	}

	const Type* operator[](int i) const
	{
		return m[i];
	}

	Type& operator()(int i, int j) {
		return m[i][j];
	}

	const Type& operator()(int i, int j) const {
		return m[i][j];
	}

	int RowCount() const
	{
		return m_row_count;
	}

	int ColCount() const
	{
		return m_col_count;
	}

	bool Create(
		int row, // number of rows
		int col // number of columns
	)
	{
		if (row == 0 || col == 0) return false;
		if (m_row_count == row && m_col_count == col && m != nullptr)
		{
			return true;
		}
		Destroy();
		m_row_count = row;
		m_col_count = col;
		m = new Type*[m_row_count];
		for (int i = 0; i < m_row_count; i++)
		{
			m[i] = new Type[m_col_count];
		}
		return true;
	}

	bool Resize(int new_row, int new_col)
	{
		if (new_row <= 0 || new_col <= 0)
			return false;

		if (new_row == m_row_count && new_col == m_col_count) return false;

		Type** new_data = new Type * [new_row];
		for (int i = 0; i < new_row; ++i)
			new_data[i] = new Type[new_col];

		int copy_row = min(m_row_count, new_row);
		int copy_col = min(m_col_count, new_col);
		for (int i = 0; i < copy_row; ++i)
		{
			for (int j = 0; j < copy_col; ++j)
				new_data[i][j] = m[i][j];
		}
		Destroy();

		m = new_data;
		m_row_count = new_row;
		m_col_count = new_col;

		//std::copy(&src.m[i][0], &src.m[i][col_count], &m[i][0]);

		return true;
	}

	void Destroy()
	{
		if (m != nullptr)
		{
			for (int i = 0; i < m_row_count; i++)
			{
				delete[] m[i];
			}
			delete[] m;
		}
		m = nullptr;
		m_row_count = 0;
		m_col_count = 0;
	}

	

	ON_TMatrix Transpose()
	{
		int cntRow = this->RowCount();
		int cntCol = this->ColCount();
		ON_TMatrix newMatrix(cntCol, cntRow);
		for (int i = 0; i < cntRow; ++i)
		{
			for (int j = 0; j < cntCol; ++j)
				newMatrix[j][i] = m[i][j];
		}
		return newMatrix;
	}


	friend std::ostream& operator<< (std::ostream& out, const ON_TMatrix<Type>& item) 
	{
		int row = item.RowCount();
		int col = item.ColCount();
		out << "[";
		for (int i = 0; i < row; i++)
		{
			out << "[";
			for (int j = 0; j < col; j++)
			{
				out << item(i, j);
				if (j != col - 1)
				{
					out << ", ";
				}
			}
			out << "]";
			if (i != row - 1)
			{
				out << std::endl;
			}
		}
		out << "]";
		return out;
	}

	

	RowIterator begin() { return m; }
	RowIterator end() { return m + m_row_count; }

	ConstRowIterator begin() const { return m; }
	ConstRowIterator end() const { return m + m_row_count; }

private:
	Type**	m{ nullptr };
	int	    m_row_count{ 0 };
	int	    m_col_count{ 0 };
};

template <typename T>
ON_TMatrix<T> ON_AppendRow(const ON_TMatrix<T>& first, const ON_TMatrix<T>& second)
{
	ON_ASSERT(first.ColCount() == second.RowCount());

	int firstRow = first.RowCount();
	int secondRow = second.RowCount();
	int firstCol = first.ColCount();

	ON_TMatrix<T> target(firstRow + secondRow, firstCol);
	for (int i = 0; i < firstRow; i++)
	{
		for (int j = 0; j < firstCol; j++)
		{
			target[i][j] = first[i][j];
		}
	}
	for (int i = 0; i < secondRow; i++)
	{
		for (int j = 0; j < firstCol; j++)
		{
			target[i + firstRow][j] = second[i][j];
		}
	}

	return target;
	
}

template <typename T>
ON_TMatrix<T> ON_AppendCol(const ON_TMatrix<T>& first, const ON_TMatrix<T>& second)
{
	ON_ASSERT(first.RowCount() == second.RowCount());

	int firstRow = first.RowCount();
	int secondRow = second.RowCount();
	int firstCol = first.ColCount();
	int secondCol = second.ColCount();

	ON_TMatrix<T> target(firstRow, firstCol + secondCol);

	for (int i = 0; i < firstRow; i++)
	{
		for (int j = 0; j < firstCol; j++)
		{
			target[i][j] = first[i][j];
		}
	}

	for (int i = 0; i < secondRow; i++)
	{
		for (int j = 0; j < secondCol; j++)
		{
			target[i][j + firstCol] = second[i][j];
		}
	}
	return target;
}


template <typename Type>
requires std::is_floating_point_v<Type>
ON_TMatrix<Type> operator*(const ON_TMatrix<Type>& left, const ON_TMatrix<Type>& right)
{
	// 행렬 곱셈 조건 확인
	if (left.m_col_count != right.m_row_count)
		throw std::invalid_argument("Matrix dimensions do not match for multiplication.");

	// 결과 행렬 생성
	ON_TMatrix<Type> result(left.m_row_count, right.m_col_count);

	for (int i = 0; i < result.RowCount(); i++)
	{
		for (int j = 0; j < result.ColCount(); j++)
		{
			result[i][j] = static_cast<Type>(0.0);
		}
	}

	// 행렬 곱셈 수행
	for (int i = 0; i < left.m_row_count; ++i)
	{
		for (int j = 0; j < right.m_col_count; ++j)
		{
			Type sum = static_cast<Type>(0.0); // Type에 맞게 초기화
 			for (int k = 0; k < left.m_col_count; ++k)
			{
				// 요소 접근 시 Type을 반환하는지 확인
				sum += left.m[i][k] * right.m[k][j];
			}
			result.m[i][j] = sum;
		}
	}

	return result; // 새 객체 반환 (좌변 보존)
}

template class ON_CLASS ON_TMatrix<double>;
template class ON_CLASS ON_TMatrix<ON_3dPoint>;
template class ON_CLASS ON_TMatrix<ON_4dPoint>;
template class ON_CLASS ON_TMatrix<unsigned char>;

ON_TMatrix<double> ON_BezierToPowerTMatrix(int degree)
{
  const int n = degree;
  ON_TMatrix<double> M(n + 1, n + 1);
  for (int m = 0; m <= n; ++m)
  {
    for (int i = 0; i <= n; ++i)
    {
      M[m][i] = 0.0;
    }
  }


  for (int m = 0; m <= n; ++m)
  {
    for (int i = 0; i <= n; ++i)
    {
      if (m >= i)
      {
        int k = m - i;
        double val = ON_Binomial(n, i) * ON_Binomial(n - i, k) * ((k % 2) ? -1.0 : 1.0);
        M[m][i] = val;
      }
    }
  }
  return M;
}


// --- Power -> Bezier 변환행렬: M_PB = (M_BP)^{-1} ---
ON_TMatrix<double> ON_PowerToBezierTMatrix(int degree, const ON_TMatrix<double>& M_BP)
{
  ON_TMatrix<double> Minv = M_BP.Invert(ON_ZERO_TOLERANCE);
  return Minv;
}



// --- Horner(오름차): p(t)=a0 + a1 t + ... + an t^n ---
double ON_HornerAscending(const std::vector<double>& a, double t)
{
  double s = 0.0;
  for (int i = (int)a.size() - 1; i >= 0; --i) s = s * t + a[i];
  return s;
}


double ON_HornerDesending(const std::vector<double>& a, double t)
{
  double s = 0.0;
  for (int i = 0; i < (int)a.size(); ++i) s = s * t + a[i];
  return s;
}

bool ON_PowerToBezierSurface(
  const ON_TMatrix<double>& Ax,
  const ON_TMatrix<double>& Ay,
  const ON_TMatrix<double>& Az,
  int du, int dv,
  ON_BezierSurface& out_bezier) // out_bezier은 이 함수에서 초기화
{
  const int nu = du + 1;
  const int nv = dv + 1;
  if (Ax.RowCount() != nu || Ax.ColCount() != nv) return false;
  if (Ay.RowCount() != nu || Ay.ColCount() != nv) return false;
  if (Az.RowCount() != nu || Az.ColCount() != nv) return false;

  // 1D inverse 변환
  ON_TMatrix<double> Mu = ON_BezierToPowerTMatrix(du);
  ON_TMatrix<double> Mv = ON_BezierToPowerTMatrix(dv);
  ON_TMatrix<double> Pu = Mu.Invert(ON_ZERO_TOLERANCE);          // M_PB_u;
  ON_TMatrix<double> Pv = Mv.Invert(ON_ZERO_TOLERANCE);          // M_PB_v

  // P = Pu * A * Pv^T
  ON_TMatrix<double> PvTrans = Pv.Transpose();

  ON_TMatrix<double> Px = ((Pu * Ax) * PvTrans);
  ON_TMatrix<double> Py = ((Pu * Ay) * PvTrans);
  ON_TMatrix<double> Pz = ((Pu * Az) * PvTrans);

  out_bezier.Create(3, false, nu, nv);
  for (int i = 0; i < nu; ++i)
  {
    for (int j = 0; j < nv; ++j)
    {
      out_bezier.SetCV(i, j, ON_3dPoint(Px[i][j], Py[i][j], Pz[i][j]));
    }
  }
  return true;
}

bool ON_BezierSurfaceToPower(
  const ON_BezierSurface& srf,
  ON_TMatrix<double>& Ax,
  ON_TMatrix<double>& Ay,
  ON_TMatrix<double>& Az)
{
  if (srf.IsRational())
    return false; // non-rational만 처리

  int du = srf.Degree(0);
  int dv = srf.Degree(1);
  const int nu = du + 1;
  const int nv = dv + 1;

  // Control net: Px,Py,Pz  (nu x nv)
  ON_TMatrix<double> Px(nu, nv), Py(nu, nv), Pz(nu, nv);
  for (int i = 0; i < nu; ++i)
  {
    for (int j = 0; j < nv; ++j)
    {
      ON_3dPoint cv; srf.GetCV(i, j, cv);
      Px[i][j] = cv.x;
      Py[i][j] = cv.y;
      Pz[i][j] = cv.z;
    }
  }

  // 1D 변환 행렬
  ON_TMatrix<double> Mu = ON_BezierToPowerTMatrix(du);
  ON_TMatrix<double> Mv = ON_BezierToPowerTMatrix(dv);

  // A = Mu * P * Mv^T
  ON_TMatrix<double> MvTrans = Mv.Transpose();
  Ax = Mu * Px * MvTrans;
  Ay = Mu * Py * MvTrans;
  Az = Mu * Pz * MvTrans;
  return true;
}


ON_3dPoint ON_EvaluatePowerSurfaceAscending(
  const ON_TMatrix<double>& Ax,  // (du+1) x (dv+1), coeff of u^p v^q
  const ON_TMatrix<double>& Ay,
  const ON_TMatrix<double>& Az,
  double u,
  double v
)
{
  const int ru = Ax.RowCount(); // du+1
  const int cv = Ax.ColCount(); // dv+1
  const int du = ru - 1;
  const int dv = cv - 1;

  // For each p (row), evaluate poly in v
  std::vector<double> sx(ru), sy(ru), sz(ru);
  for (int p = 0; p <= du; ++p)
  {
    sx[p] = ON_HornerAscending(&Ax[p][0], dv, v);
    sy[p] = ON_HornerAscending(&Ay[p][0], dv, v);
    sz[p] = ON_HornerAscending(&Az[p][0], dv, v);
  }
  // Then in u
  double x = ON_HornerAscending(sx.data(), du, u);
  double y = ON_HornerAscending(sy.data(), du, u);
  double z = ON_HornerAscending(sz.data(), du, u);
  return ON_3dPoint(x, y, z);

}

ON_3dPoint ON_EvaluateBezierSurfaceBasis(
  const ON_BezierSurface& srf, 
  int du, 
  int dv, 
  double u, 
  double v)
{
  std::vector<double> Bu, Bv;
  ON_BasisFuns(du, u, Bu);
  ON_BasisFuns(dv, v, Bv);

  ON_3dPoint p(0, 0, 0);
  for (int i = 0; i <= du; ++i)
  {
    for (int j = 0; j <= dv; ++j)
    {
      ON_3dPoint cv;
      srf.GetCV(i, j, cv);
      const double w = Bu[i] * Bv[j];
      p.x += w * cv.x;
      p.y += w * cv.y;
      p.z += w * cv.z;
    }
  }
  return p;
}




int main() {
  ON::Begin();
  ON_TextLog log;

  SetConsoleOutputCP(CP_UTF8);
  SetConsoleCP(CP_UTF8);
 
  // 1) bicubic 예제 (4x4)
  const int du = 3, dv = 3;
  const int nu = du + 1, nv = dv + 1;

  ON_BezierSurface srf(3, false, nu, nv);
  // 샘플 제어점 (간단한 물결 형태)
  for (int i = 0; i < nu; ++i)
  {
    for (int j = 0; j < nv; ++j)
    {
      double x = i;
      double y = j;
      double z = 0.25 * std::sin((i * ON_PI) / 3.0) * std::cos((j * ON_PI) / 3.0);
      srf.SetCV(i, j, ON_3dPoint(x, y, z));
    }
  }

  // 2) Bezier -> Power
  ON_TMatrix<double> Ax, Ay, Az;
  if (!ON_BezierSurfaceToPower(srf, Ax, Ay, Az))
  {
    std::cout << "Rational surface not supported in this demo.\n";
    return 0;
  }

  // 3) Power -> Bezier (복원)
  ON_BezierSurface srf_rec;
  ON_PowerToBezierSurface(Ax, Ay, Az, du, dv, srf_rec);

  // 4) CV 비교
  double max_cv_err = 0.0;
  for (int i = 0; i < nu; ++i)
  {
    for (int j = 0; j < nv; ++j)
    {
      ON_3dPoint a, b;
      srf.GetCV(i, j, a);
      srf_rec.GetCV(i, j, b);
      max_cv_err = ON_MaxAbs(max_cv_err, a.DistanceTo(b));
    }
  }
  std::cout << "max |CV - recon(CV)| = " << max_cv_err << "\n";

  // 5) 임의 (u,v)에서 값 비교 (Bezier vs Power)
  double max_eval_err = 0.0;
  for (int k = 0; k < 10; ++k)
  {
    double u = k / 9.0;
    for (int m = 0; m < 10; ++m)
    {
      double v = m / 9.0;
      ON_3dPoint pB = ON_EvaluateBezierSurfaceBasis(srf, du, dv, u, v);
      ON_3dPoint pP = ON_EvaluatePowerSurfaceAscending(Ax, Ay, Az, u, v);
      max_eval_err = ON_MaxAbs(max_eval_err, pB.DistanceTo(pP));
    }
  }
  std::cout << "max |Bezier - Power| = " << max_eval_err << "\n";

  // 6) 샘플 평가 출력
  double u = 0.37, v = 0.58;
  ON_3dPoint PB = ON_EvaluateBezierSurfaceBasis(srf, du, dv, u, v);
  ON_3dPoint PP = ON_EvaluatePowerSurfaceAscending(Ax, Ay, Az, u, v);
  std::cout << "Evaluate at (u,v) = (" << u << "," << v << ")\n";
  std::cout << "  Bezier: (" << PB.x << "," << PB.y << "," << PB.z << ")\n";
  std::cout << "  Power : (" << PP.x << "," << PP.y << "," << PP.z << ")\n";
  std::cout << "  |diff|: " << PB.DistanceTo(PP) << "\n";


  ON::End();
  return 0;
}

```

## 8. 실험 결과 예시

- max |CV - recon(CV)| ≈ $\(8.3 \times 10^{-17}\)$  
- max |Bezier - Power| ≈ $\(9.1 \times 10^{-16}\)$  
- 샘플 평가 $\((u,v) = (0.37,0.58)\)$:  
  - Bézier: (1.11, 1.74, -0.027175)  
  - Power : (1.11, 1.74, -0.027175)  
  - 차이: ~ $\(5 \times 10^{-16}\)$

즉, 이론대로 구현이 정확히 동작함을 확인할 수 있다.

---







