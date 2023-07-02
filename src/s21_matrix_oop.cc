// created by pizpotli
#include "s21_matrix_oop.h"

// KONSTRUCTORS

S21Matrix::S21Matrix() noexcept : rows_(0), cols_(0) { matrix_ = nullptr; }

S21Matrix::S21Matrix(int rows, int cols) {
  MallocMatrix(rows, cols);
  ZeroMatrix();
}

S21Matrix::S21Matrix(const S21Matrix& other) {
  MallocMatrix(other.rows_, other.cols_);
  CopyMatrix(other);
}

S21Matrix::S21Matrix(S21Matrix&& other) noexcept
    : rows_(other.rows_), cols_(other.cols_) {
  if (this != &other) {
    matrix_ = other.matrix_;
    other.matrix_ = nullptr;
    other.rows_ = 0;
    other.cols_ = 0;
  }
}

// DESTRUCTOR

S21Matrix::~S21Matrix() noexcept { Remove(); }

// ARITHMETICS

void S21Matrix::SumMatrix(const S21Matrix& other) {
  CheckMistakes(other, 1);
  CheckMistakes(other, 2);
  for (int i = 0; i < other.rows_; i++) {
    for (int j = 0; j < other.cols_; j++) {
      matrix_[i][j] += other.matrix_[i][j];
    }
  }
}

bool S21Matrix::EqMatrix(const S21Matrix& other) const noexcept {
  bool result = true;
  if (rows_ == other.rows_ && cols_ == other.cols_) {
    for (int i = 0; i < rows_ && result == true; i++) {
      for (int j = 0; j < cols_ && result == true; j++) {
        if (fabs(matrix_[i][j] - other.matrix_[i][j]) > 1e-6) {
          result = false;
        }
      }
    }
  } else {
    result = false;
  }
  return result;
}

void S21Matrix::SubMatrix(const S21Matrix& other) {
  CheckMistakes(other, 1);
  CheckMistakes(other, 2);
  for (int i = 0; i < other.rows_; i++) {
    for (int j = 0; j < other.cols_; j++) {
      matrix_[i][j] -= other.matrix_[i][j];
    }
  }
}

void S21Matrix::MulNumber(const double num) {
  CheckMistakes2(1);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] *= num;
    }
  }
}

void S21Matrix::MulMatrix(const S21Matrix& other) {
  CheckMistakes(other, 1);
  CheckMistakes(other, 3);
  S21Matrix tmp(rows_, cols_);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < other.cols_; j++) {
      double res = 0;
      for (int k = 0; k < cols_; k++) {
        res += matrix_[i][k] * other.matrix_[k][j];
      }
      tmp.matrix_[i][j] = res;
    }
  }
  *this = tmp;
  tmp.Remove();
}

S21Matrix S21Matrix::Transpose() const {
  CheckMistakes2(1);
  S21Matrix tmp(cols_, rows_);
  for (int i = 0; i < tmp.rows_; i++) {
    for (int j = 0; j < tmp.cols_; j++) {
      tmp(i, j) = matrix_[j][i];
    }
  }
  return tmp;
}

S21Matrix S21Matrix::CalcComplements() const {
  CheckMistakes2(1);
  CheckMistakes2(2);
  S21Matrix result(*this);
  S21Matrix minor(*this);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      S21Matrix tmp(rows_ - 1, cols_ - 1);
      minor.Minor(i, j, tmp);
      double det = 0;
      det = tmp.Determinant();
      result.matrix_[i][j] = pow(-1.0, i + j) * det;
      tmp.Remove();
    }
  }
  minor.Remove();
  return result;
}

S21Matrix S21Matrix::InverseMatrix() const {
  CheckMistakes2(1);
  double det = 0;
  S21Matrix tmp(*this);
  det = this->Determinant();
  if (fabs(det) < 1e-7) {
    tmp.Remove();
    throw std::out_of_range("ERROR: calculation impossible: Determinant = 0");
  } else {
    if (rows_ > 1 && cols_ > 1) {
      S21Matrix tmp2 = tmp.CalcComplements();
      S21Matrix transpon = tmp2.Transpose();
      transpon.MulNumber(1.0 / det);
      tmp.CopyMatrix(transpon);
      transpon.Remove();
      tmp2.Remove();
    } else {
      S21Matrix tmp2(1, 1);
      tmp.matrix_[0][0] = 1.0 / matrix_[0][0];
    }
  }
  return tmp;
}

double S21Matrix::Determinant() const {
  CheckMistakes2(1);
  CheckMistakes2(2);
  double res = 0;
  S21Matrix tmp;
  S21Matrix copy(*this);
  res = (-copy.Triangulate(tmp) % 2) ? -1 : 1;
  for (int x = 0; x < tmp.cols_; x++) {
    res *= tmp.matrix_[x][x];
  }
  tmp.Remove();
  copy.Remove();
  return -res;
}

// OPERATORS

S21Matrix S21Matrix::operator+(const S21Matrix& other) const {
  S21Matrix res(*this);
  res.SumMatrix(other);
  return res;
}

S21Matrix S21Matrix::operator-(const S21Matrix& other) const {
  S21Matrix res(*this);
  res.SubMatrix(other);
  return res;
}

S21Matrix S21Matrix::operator*(const S21Matrix& other) const {
  S21Matrix res(*this);
  res.MulMatrix(other);
  return res;
}

S21Matrix S21Matrix::operator*(const double number) const {
  S21Matrix res(*this);
  res.MulNumber(number);
  return res;
}

S21Matrix operator*(const double number, const S21Matrix& other) {
  S21Matrix tmp(other);
  return tmp *= number;
}

S21Matrix& S21Matrix::operator+=(const S21Matrix& other) {
  SumMatrix(other);
  return *this;
}

S21Matrix& S21Matrix::operator-=(const S21Matrix& other) {
  SubMatrix(other);
  return *this;
}

S21Matrix& S21Matrix::operator*=(const S21Matrix& other) {
  MulMatrix(other);
  return *this;
}

S21Matrix& S21Matrix::operator*=(const double number) {
  MulNumber(number);
  return *this;
}

bool S21Matrix::operator==(const S21Matrix& other) const {
  bool res = EqMatrix(other);
  return res;
}

S21Matrix& S21Matrix::operator=(const S21Matrix& other) {
  CopyMatrix(other);
  return *this;
}

double& S21Matrix::operator()(const int x, const int y) {
  if (x >= rows_ || y >= cols_ || x < 0 || y < 0) {
    throw std::out_of_range("ERROR: index outside matrix");
  }
  return matrix_[x][y];
}

// ACCESSORS

int S21Matrix::GetRows() const noexcept { return rows_; }

int S21Matrix::GetCols() const noexcept { return cols_; }

// MUTATORS

void S21Matrix::SetRows(const int rows) {
  if (rows < 1) {
    throw std::out_of_range("ERROR: incorrect matrix");
  }
  S21Matrix A;
  A.MallocMatrix(rows, cols_);
  A.ZeroMatrix();
  int a = 0;
  if (A.rows_ < rows_) {
    a = A.rows_;
  } else {
    a = rows_;
  }
  for (int i = 0; i < a; i++) {
    for (int j = 0; j < A.cols_; j++) {
      A.matrix_[i][j] = matrix_[i][j];
    }
  }
  CopyMatrix(A);
  A.Remove();
}

void S21Matrix::SetCols(const int cols) {
  if (cols < 1) {
    throw std::out_of_range("ERROR: incorrect matrix");
  }
  S21Matrix A;
  A.MallocMatrix(rows_, cols);
  A.ZeroMatrix();
  int a = 0;
  if (A.cols_ < cols_) {
    a = A.cols_;
  } else {
    a = cols_;
  }
  for (int i = 0; i < A.rows_; i++) {
    for (int j = 0; j < a; j++) {
      A.matrix_[i][j] = matrix_[i][j];
    }
  }
  CopyMatrix(A);
  A.Remove();
}

// HELP FUNCTIONS

void S21Matrix::Remove() noexcept {
  if (matrix_) {
    for (int i = 0; i < rows_; i++) {
      if (matrix_[i]) {
        delete [] matrix_[i];
      }
    }
    delete [] matrix_;
    matrix_ = nullptr;
  }
  rows_ = 0;
  cols_ = 0;
}

void S21Matrix::MallocMatrix(int x, int y) {
  if (x < 1 || y < 1) {
    throw std::out_of_range("ERROR: incorrect matrix");
  }
  matrix_ = new double*[x]();
  for (int i = 0; i < x; i++) {
    matrix_[i] = new double[y]();
  }
  rows_ = x;
  cols_ = y;
}

void S21Matrix::SumStr(int row, int ro, double tmp) noexcept {
  for (int x = 0; x < cols_; x++) matrix_[row][x] += matrix_[ro][x] * tmp;
}

void S21Matrix::Minor(int x, int y, S21Matrix& other) noexcept {
  int i1 = 0, i2 = 0, j1 = 0, j2 = 0;
  for (i1 = 0; i1 < rows_ - 1; i1++) {
    if (i1 == x) {
      i2 = 1;
    }
    j2 = 0;
    for (j1 = 0; j1 < rows_ - 1; j1++) {
      if (j1 == y) {
        j2 = 1;
      }
      other.matrix_[i1][j1] = matrix_[i1 + i2][j1 + j2];
    }
  }
}

void S21Matrix::CopyMatrix(const S21Matrix& other) noexcept {
  for (int i = 0; i < rows_; i++) {
    delete (matrix_[i]);
  }
  delete (matrix_);
  matrix_ = new double*[other.rows_]();
  for (int i = 0; i < other.rows_; i++) {
    matrix_[i] = new double[other.cols_]();
  }
  rows_ = other.rows_;
  cols_ = other.cols_;
  for (int i = 0; i < other.rows_; i++) {
    for (int j = 0; j < other.cols_; j++) {
      matrix_[i][j] = other.matrix_[i][j];
    }
  }
}

void S21Matrix::ZeroMatrix() noexcept {
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = 0;
    }
  }
}

int S21Matrix::Sdvig(int i, int j) noexcept {
  double* tmp = matrix_[i];
  matrix_[i] = matrix_[j];
  matrix_[j] = tmp;
  return -1;
}

int S21Matrix::Triangulate(S21Matrix& other) noexcept {
  int znak = 1;
  other.CopyMatrix(*this);
  for (int i = 0; i < other.cols_; i++) {
    if (!other.matrix_[i][i]) {
      for (int j = i + 1; j < other.rows_; j++) {
        if (other.matrix_[j][i]) {
          znak -= other.Sdvig(j, i);
        }
      }
    }
    for (int j = i + 1; j < rows_; j++) {
      other.SumStr(j, i, -other.matrix_[j][i] / other.matrix_[i][i]);
    }
  }
  return -znak;
}

void S21Matrix::CheckMistakes(const S21Matrix& other, const int number) const {
  if (number == 1) {
    if (rows_ < 1 || cols_ < 1 || other.rows_ < 1 || other.cols_ < 1) {
      throw std::out_of_range("ERROR: incorrect matrix");
    }
  } else if (number == 2) {
    if (rows_ != other.rows_ || cols_ != other.cols_) {
      throw std::out_of_range("ERROR: different dimensions of matrices");
    }
  } else if (number == 3) {
    if (cols_ != other.rows_) {
      throw std::out_of_range("ERROR: sides are not equal");
    }
  }
}

void S21Matrix::CheckMistakes2(const int number) const {
  if (number == 1) {
    if (rows_ < 1 || cols_ < 1) {
      throw std::out_of_range("ERROR: incorrect matrix");
    }
  } else if (number == 2) {
    if (cols_ != rows_) {
      throw std::out_of_range("ERROR: matrix is not square");
    }
  }
}
