// created by pizpotli
#ifndef CPP1_S21_MATRIXPLUS_3_SRC_S21_MATRIX_OOP_H_
#define CPP1_S21_MATRIXPLUS_3_SRC_S21_MATRIX_OOP_H_

#include <cmath>
#include <iostream>

class S21Matrix {
 public:
  // Konstructors

  S21Matrix() noexcept;
  S21Matrix(int rows, int cols);
  S21Matrix(const S21Matrix& other);
  S21Matrix(S21Matrix&& other) noexcept;

  // Destructor

  ~S21Matrix() noexcept;

  // Arithmetics

  void SumMatrix(const S21Matrix& other);
  bool EqMatrix(const S21Matrix& other) const noexcept;
  void SubMatrix(const S21Matrix& other);
  void MulNumber(const double num);
  void MulMatrix(const S21Matrix& other);
  S21Matrix Transpose() const;
  S21Matrix CalcComplements() const;
  S21Matrix InverseMatrix() const;
  double Determinant() const;

  // operators

  S21Matrix operator+(const S21Matrix& other) const;
  S21Matrix operator-(const S21Matrix& other) const;
  S21Matrix operator*(const S21Matrix& other) const;
  S21Matrix operator*(const double number) const;
  bool operator==(const S21Matrix& other) const;

  S21Matrix& operator+=(const S21Matrix& other);
  S21Matrix& operator-=(const S21Matrix& other);
  S21Matrix& operator*=(const S21Matrix& other);
  S21Matrix& operator*=(const double number);
  S21Matrix& operator=(const S21Matrix& other);

  double& operator()(const int x, const int y);
  friend S21Matrix operator*(const double number, const S21Matrix& other);

  // Accessors

  int GetRows() const noexcept;
  int GetCols() const noexcept;

  // Mutators

  void SetRows(const int rows);
  void SetCols(const int cols);

 private:
  int rows_, cols_;
  double** matrix_;

  // help functions

  void Remove() noexcept;
  void MallocMatrix(int x, int y);
  void SumStr(int row, int ro, double tmp) noexcept;
  void Minor(int i, int j, S21Matrix& other) noexcept;
  void CopyMatrix(const S21Matrix& other) noexcept;
  void ZeroMatrix() noexcept;
  int Sdvig(int i, int j) noexcept;
  int Triangulate(S21Matrix& other) noexcept;
  void CheckMistakes(const S21Matrix& other, const int number) const;
  void CheckMistakes2(const int number) const;
};

S21Matrix operator*(const double number, const S21Matrix& other);

#endif  // CPP1_S21_MATRIXPLUS_3_SRC_S21_MATRIX_OOP_H_