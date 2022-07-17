//
//  NumMethods.c
//
//
//  Created by Andrea Vivani on 11/08/16.
//  Copyright © 2016 Andrea Vivani. All rights reserved.
//

#include "NumMethods.h"
#include "math.h"

//-------------------Forward substitution----------------------//
// assumes that the matrix A is already a lower triangular one. No check!

void fwsub(Matrix A, Matrix B, Matrix result)
{
  int16_t i, j, k;
  float tmp;
  for (k = 0; k < B.cols; k++) {
    ELEM(result, 0, k) = ELEM(B, 0, k) / ELEM(A, 0, 0);
    for (i = 1; i < A.rows; i++) {
      tmp = 0.0;
      for (j = 0; j < i; j++) {
        tmp += ELEM(A, i, j) * (ELEM(result, j, k));
      }
      ELEM(result, i, k) = (ELEM(B, i, k) - tmp) / ELEM(A, i, i);
    }
  }
  return;
}

//---------------Forward substitution with permutation-------------------//
// assumes that the matrix A is already a lower triangular one. No check!

void fwsubPerm(Matrix A, Matrix B, Matrix P, Matrix result)
{
  int16_t i, j, k;
  float tmp;
  for (k = 0; k < B.cols; k++) {
    ELEM(result, 0, k) = ELEM(B, (uint8_t) ELEM(P, 0, 0), k) / ELEM(A, 0, 0);
    for (i = 1; i < A.rows; i++) {
      tmp = 0.0;
      for (j = 0; j < i; j++) {
        tmp += ELEM(A, i, j) * ELEM(result, j, k);
      }
      ELEM(result, i, k) = (ELEM(B, (uint8_t) ELEM(P, i, 0), k) - tmp) / ELEM(A, i, i);
    }
  }
  return;
}

//-------------------Backward substitution----------------------//
// assumes that the matrix A is already an upper triangular one. No check!

void bksub(Matrix A, Matrix B, Matrix result)
{
  int16_t i, j, k;
  float tmp;
  for (k = 0; k < B.cols; k++) {
    ELEM(result, A.cols - 1, k) = ELEM(B, A.cols - 1, k) / ELEM(A, A.cols - 1, A.cols - 1);
    for (i = A.rows - 2; i >= 0; i--) {
      tmp = 0.0;
      for (j = A.cols - 1; j > i; j--) {
        tmp += ELEM(A, i, j) * ELEM(result, j, k);
      }
      ELEM(result, i, k) = (ELEM(B, i, k) - tmp) / ELEM(A, i, i);
    }
  }
  return;
}

//--------------Backward substitution with permutation-----------------//
// assumes that the matrix A is already an upper triangular one. No check!

void bksubPerm(Matrix A, Matrix B, Matrix P, Matrix result)
{
  int16_t i, j, k;
  float tmp;
  for (k = 0; k < B.cols; k++) {
    ELEM(result, A.cols - 1, k) = ELEM(B, (uint8_t) ELEM(P, A.cols - 1, 0), k) / ELEM(A, A.cols - 1, A.cols - 1);
    for (i = A.rows - 2; i >= 0; i--) {
      tmp = 0.0;
      for (j = A.cols - 1; j > i; j--) {
        tmp += ELEM(A, i, j) * ELEM(result, j, k);
      }
      ELEM(result, i, k) = (ELEM(B, (uint8_t) ELEM(P, i, 0), k) - tmp) / ELEM(A, i, i);
    }
  }
  return;
}

//-------------------------LU factorization using Crout's Method--------------------------------//
// factorizes the A matrix as the product of a unit upper triangular matrix U and a lower triangular matrix L

uint8_t LU_Crout(Matrix A, Matrix L, Matrix U)
{
  int16_t ii, jj, kk;
  float sum = 0.0;
  matIdentity(U);
  matZeros(L);
  for (jj = 0; jj < A.rows; jj++) {
    for (ii = jj; ii < A.rows; ii++) {
      sum = 0.0f;
      for (kk = 0; kk < jj; kk++) {
        sum += ELEM(L, ii, kk) * ELEM(U, kk, jj);
      }
      ELEM(L, ii, jj) = ELEM(A, ii, jj) - sum;
    }

    for (ii = jj; ii < A.rows; ii++) {
      sum = 0;
      for (kk = 0; kk < jj; kk++) {
        sum += ELEM(L, jj, kk) * ELEM(U, kk, ii);
      }
      if (ELEM(L, jj, jj) == 0) {
        return 0;
      }
      ELEM(U, jj, ii) = (ELEM(A, jj, ii) - sum) / ELEM(L, jj, jj);
    }
  }
  return 1;
}

//-------------------------LU factorization using Cormen's Method--------------------------------//
// factorizes the A matrix as the product of a unit upper triangular matrix U and a lower triangular matrix L

uint8_t LU_Cormen(Matrix A, Matrix L, Matrix U)
{
  int16_t i, j, k;
  float tmp;
  float _A_cp_data[A.rows * A.cols];
  Matrix A_cp = newMatrix(A.rows, A.cols, _A_cp_data);
  matCopyData(A, A_cp);
  matZeros(U);
  matIdentity(L);

  for (k = 0; k < A_cp.rows; k++) {
    ELEM(U, k, k) = ELEM(A_cp, k, k);
    if (ELEM(A_cp, k,k) == 0) {
      return 0;
    }
    tmp = 1.0 / ELEM(U, k, k);
    for (i = k + 1; i < A_cp.rows; i++) {
      ELEM(L, i, k) = ELEM(A_cp, i, k) * tmp;
      ELEM(U, k, i) = ELEM(A_cp, k, i);
    }
    for (i = k + 1; i < A_cp.rows; i++) {
      for (j = k + 1; j < A_cp.rows; j++) {
        ELEM(A_cp, i, j) -= ELEM(L, i, k) * ELEM(U, k, j);
      }
    }
  }
  return 1;
}

//-----------------------LUP factorization using Cormen's Method------------------------------//
// factorizes the A matrix as the product of a upper triangular matrix U and a unit lower triangular matrix L
// returns the factor that has to be multiplied to the determinant of U in order to obtain the correct value

int8_t LUP_Cormen(Matrix A, Matrix L, Matrix U, Matrix P)
{
  int16_t i, j, k;
  float tmp, tmp2;
  int16_t pivrow;
  int8_t d_mult = 1; // determinant multiplying factor
  float _A_cp_data[A.rows * A.cols];
  Matrix A_cp = newMatrix(A.rows, A.cols, _A_cp_data);
  matCopyData(A, A_cp);
  matZeros(U);
  matIdentity(L);
  // initialization
  for (i = 0; i < A_cp.rows; i++) {
    ELEM(P, i, 0) = i;
  }

  // outer loop over diagonal pivots
  for (k = 0; k < A_cp.rows - 1; k++) {

    // inner loop to find the largest pivot
    pivrow = k;
    tmp = fabs(ELEM(A_cp, k, k));
    for (i = k + 1; i < A_cp.rows; i++) {
      tmp2 = fabs(ELEM(A_cp, i, k));
      if (tmp2 > tmp) {
        tmp = tmp2;
        pivrow = i;
      }
    }
    // check for singularity
    if (ELEM(A_cp, pivrow, k) == 0) {
      return 0;
    }

    // swap rows
    if (pivrow != k) {
      tmp = ELEM(P, k, 0);
      ELEM(P, k, 0) = ELEM(P, pivrow, 0);
      ELEM(P, pivrow, 0) = tmp;
      d_mult *= -1;

      for (j = 0; j < A_cp.rows; j++) {
        tmp = ELEM(A_cp, k, j);
        ELEM(A_cp, k, j) = ELEM(A_cp, pivrow, j);
        ELEM(A_cp, pivrow, j) = tmp;
      }
    }
    tmp = 1.0 / ELEM(A_cp, k, k);
    // Gaussian elimination
    for (i = k + 1; i < A_cp.rows; i++) { // iterate down rows
      ELEM(A_cp, i, k) *= tmp;
      for (j = k + 1; j < A_cp.rows; j++) { // iterate across rows
        ELEM(A_cp, i, j) -= ELEM(A_cp, i, k) * ELEM(A_cp, k, j);
      }
    }
  }
  for (k = 0; k < A_cp.rows; k++) {
    ELEM(U, k, k) = ELEM(A_cp, k, k);
    for (j = k + 1; j < A_cp.rows; j++) {
      ELEM(L, j, k) = ELEM(A_cp, j, k);
      ELEM(U, k, j) = ELEM(A_cp, k, j);
    }
  }
  return d_mult;
}

//-----------------------Linear system solver using LU factorization---------------------------//
// solves the linear system A*X=B, where A is a n-by-n matrix and B an n-by-m matrix, giving the n-by-m matrix X

void LinSolveLU(Matrix A, Matrix B, Matrix result)
{
  float _L_data[A.rows * A.cols];
  float _U_data[A.cols * A.cols];
  Matrix L = newMatrix(A.rows, A.cols, _L_data);
  Matrix U = newMatrix(A.cols, A.cols, _U_data);
  //Matrix tmp1 = newMatrix(A.rows, B.cols);
  LU_Cormen(A, L, U);
  //fwsub(L, B, tmp1);
  //bksub(U, tmp1, result);
  fwsub(L, B, result);
  bksub(U, result, result); //hope it can work in-place
  return;
}

//----------------------Linear system solver using LUP factorization--------------------------//
// solves the linear system A*X=B, where A is a n-by-n matrix and B an n-by-m matrix, giving the n-by-m matrix X

void LinSolveLUP(Matrix A, Matrix B, Matrix result)
{
  float _L_data[A.rows * A.cols];
  float _U_data[A.cols * A.cols];
  float _P_data[A.rows];
  float _tmp_data[A.rows * B.cols];
  Matrix L = newMatrix(A.rows, A.cols, _L_data);
  Matrix U = newMatrix(A.cols, A.cols, _U_data);
  Matrix P = newMatrix(A.rows, 1, _P_data);
  Matrix tmp = newMatrix(A.rows, B.cols, _tmp_data);

  LUP_Cormen(A, L, U, P);
  fwsubPerm(L, B, P, tmp);
  bksub(U, tmp, result);
  return;
}

//------------Linear system solver using Gauss elimination with partial pivoting---------------//
// solves the linear system A*X=B, where A is a n-by-n matrix and B an n-by-m matrix, giving the n-by-m matrix X

void LinSolveGauss(Matrix A, Matrix B, Matrix result)
{
  uint8_t pivrow = 0;     // keeps track of current pivot row
  uint8_t k, i, j; // k: overall index along diagonals; i: row index; j: col index
  float tmp;      // used for finding max value and making row swaps
  float tmp2; // used to store abs when finding max value and to store coefficient value when eliminating values
  float _A_cp_data[A.rows * A.cols];
  Matrix A_cp = newMatrix(A.rows, A.cols, _A_cp_data);
  matCopyData(A, A_cp);
  float _B_cp_data[B.rows * B.cols];
  Matrix B_cp = newMatrix(B.rows, B.cols, _B_cp_data);
  matCopyData(B, B_cp);

  for (k = 0; k < (A_cp.cols - 1); k++) {

    // find pivot row, the row with biggest entry in current column
    tmp = fabs(ELEM(A_cp, k, k));
    pivrow = k;
    for (i = k + 1; i < A_cp.cols; i++) {
      tmp2 = fabs(ELEM(A_cp, i, k)); // 'Avoid using other functions inside abs()?'
      if (tmp2 > tmp) {
        tmp = tmp2;
        pivrow = i;
      }
    }

    // check for singular Matrix
    if (ELEM(A_cp, pivrow, k) == 0.0) {
      matZeros(result);
      return;
    }

    // Execute pivot (row swap) if needed
    if (pivrow != k) {
      // swap row k of matrix A with pivrow
      for (j = k; j < A_cp.cols; j++) {
        tmp = ELEM(A_cp, k, j);
        ELEM(A_cp, k, j) = ELEM(A_cp, pivrow, j);
        ELEM(A_cp, pivrow, j) = tmp;
      }
      // swap row k of matrix B with pivrow
      for (j = 0; j < B_cp.cols; j++) {
        tmp = ELEM(B_cp, k, j);
        ELEM(B_cp, k, j) = ELEM(B_cp, pivrow, j);
        ELEM(B_cp, pivrow, j) = tmp;
      }
    }

    // Row reduction
    tmp = 1.0 / ELEM(A_cp, k, k);    // invert pivot element
    for (i = k + 1; i < A_cp.cols; i++) { // along rows
      tmp2 = ELEM(A_cp, i, k) * tmp;
      // Perform row reduction of A
      for (j = k + 1; j < A_cp.cols; j++) { //along columns of A
        ELEM(A_cp, i, j) -= tmp2 * ELEM(A_cp, k, j);
      }
      // Perform row reduction of B
      for (j = 0; j < B_cp.cols; j++) { //along columns of B
        ELEM(B_cp, i, j) -= tmp2 * ELEM(B_cp, k, j);
      }
    }

  }
  bksub(A_cp, B_cp, result);
  return;
}

//------------Gauss-Newton sensors calibration with 9 parameters---------------//
// approximates Data to a sphere of radius k by calculating 6 gains (s) and 3 biases (b), useful to calibrate some sensors (meas_sphere=S*(meas-B) with S symmetric)
// Data has n>=9 rows corresponding to the number of measures and 3 columns corresponding to the 3 axis
// X0 is the starting guess vector (usually [0 0 0 1 0 0 1 0 1]), nmax the maximum number of iterations (200 is generally fine, even if it usually converges within 10 iterations), and tol the stopping tolerance (1e-6 is usually more than fine)
/*b1=out(0,0);
 b2=out(1,0);
 b3=out(2,0);
 s11=out(3,0);
 s12=out(4,0);
 s13=out(5,0);
 s22=out(6,0);
 s23=out(7,0);
 s33=out(8,0);*/

uint8_t GaussNewton_Sens_Cal_9(Matrix Data, float k, Matrix X0, uint16_t nmax, float tol, Matrix result)
{
  matCopyData(X0, result);
  float d1, d2, d3, rx1, rx2, rx3, t1, t2, t3;
  float k2 = k * k;
  uint16_t n_iter;
  uint8_t jj;
  float _Jr_data[Data.rows * 9];
  float _res_data[Data.rows];
  float _delta_data[9];
  float _tmp1_data[9 * Data.rows];
  Matrix Jr = newMatrix(Data.rows, 9, _Jr_data);
  Matrix res = newMatrix(Data.rows, 1, _res_data);
  Matrix delta = newMatrix(9, 1, _delta_data);
  Matrix tmp1 = newMatrix(9, Data.rows, _tmp1_data);

  if ((Data.rows < 9) || (Data.cols != 3)) {
    return 0;
  }

  matZeros(Jr);
  matZeros(res);
  matZeros(delta);
  matZeros(tmp1);

  for (n_iter = 0; n_iter < nmax; n_iter++) {
    for (jj = 0; jj < Data.rows; jj++) {
      d1 = ELEM(Data, jj, 0) - ELEM(result, 0, 0);
      d2 = ELEM(Data, jj, 1) - ELEM(result, 1, 0);
      d3 = ELEM(Data, jj, 2) - ELEM(result, 2, 0);
      rx1 = -2 * (ELEM(result, 3, 0) * d1 + ELEM(result, 4, 0) * d2 + ELEM(result, 5, 0) * d3);
      rx2 = -2 * (ELEM(result, 4, 0) * d1 + ELEM(result, 6, 0) * d2 + ELEM(result, 7, 0) * d3);
      rx3 = -2 * (ELEM(result, 5, 0) * d1 + ELEM(result, 7, 0) * d2 + ELEM(result, 8, 0) * d3);
      ELEM(Jr, jj, 0) = ELEM(result, 3, 0) * rx1 + ELEM(result, 4, 0) * rx2 + ELEM(result, 5, 0) * rx3;
      ELEM(Jr, jj, 1) = ELEM(result, 4, 0) * rx1 + ELEM(result, 6, 0) * rx2 + ELEM(result, 7, 0) * rx3;
      ELEM(Jr, jj, 2) = ELEM(result, 5, 0) * rx1 + ELEM(result, 7, 0) * rx2 + ELEM(result, 8, 0) * rx3;
      ELEM(Jr, jj, 3) = -d1 * rx1;
      ELEM(Jr, jj, 4) = -d2 * rx1 - d1 * rx2;
      ELEM(Jr, jj, 5) = -d3 * rx1 - d1 * rx3;
      ELEM(Jr, jj, 6) = -d2 * rx2;
      ELEM(Jr, jj, 7) = -d3 * rx2 - d2 * rx3;
      ELEM(Jr, jj, 8) = -d3 * rx3;
      t1 = ELEM(result, 3, 0) * d1 + ELEM(result, 4,0) * d2 + ELEM(result, 5, 0) * d3;
      t2 = ELEM(result, 4, 0) * d1 + ELEM(result, 6,0) * d2 + ELEM(result, 7, 0) * d3;
      t3 = ELEM(result, 5, 0) * d1 + ELEM(result, 7,0) * d2 + ELEM(result, 8, 0) * d3;
      ELEM(res, jj, 0) = t1 * t1 + t2 * t2 + t3 * t3 - k2;
    }
    matPseudo_inv(Jr, tmp1);
    matMult(tmp1, res, delta);
    matSub(result, delta, result);
    if (matNorm(delta) < tol) {
      return 1;
    }
  }
  return 1;
}

//------------Gauss-Newton sensors calibration with 6 parameters---------------//
// approximates Data to a sphere of radius k by calculating 3 gains (s) and 3 biases (b), useful to calibrate some sensors (meas_sphere=S*(meas-B) with S diagonal)
// Data has n>=6 rows corresponding to the number of measures and 3 columns corresponding to the 3 axis
// X0 is the starting guess vector (usually [0 0 0 1 1 1]), nmax the maximum number of iterations (200 is generally fine, even if it usually converges within 10 iterations), and tol the stopping tolerance (1e-6 is usually more than fine)
/*b1=out(0,0);
 b2=out(1,0);
 b3=out(2,0);
 s11=out(3,0);
 s22=out(4,0);
 s33=out(5,0);*/

uint8_t GaussNewton_Sens_Cal_6(Matrix Data, float k, Matrix X0, uint16_t nmax, float tol, Matrix result)
{
  matCopyData(X0, result);
  float d1, d2, d3, t1, t2, t3;
  float k2 = k * k;
  uint16_t n_iter;
  uint8_t jj;
  float _Jr_data[Data.rows * 6];
  float _res_data[Data.rows];
  float _delta_data[6];
  float _tmp1_data[6 * Data.rows];
  Matrix Jr = newMatrix(Data.rows, 6, _Jr_data);
  Matrix res = newMatrix(Data.rows, 1, _res_data);
  Matrix delta = newMatrix(6, 1, _delta_data);
  Matrix tmp1 = newMatrix(6, Data.rows, _tmp1_data);

  if ((Data.rows < 6) || (Data.cols != 3)) {
    return 0;
  }

  matZeros(Jr);
  matZeros(res);
  matZeros(delta);
  matZeros(tmp1);

  for (n_iter = 0; n_iter < nmax; n_iter++) {
    for (jj = 0; jj < Data.rows; jj++) {
      d1 = ELEM(Data, jj, 0) - ELEM(result, 0, 0);
      d2 = ELEM(Data, jj, 1) - ELEM(result, 1, 0);
      d3 = ELEM(Data, jj, 2) - ELEM(result, 2, 0);
      ELEM(Jr, jj, 0) = -2 * d1 * ELEM(result, 3, 0) * ELEM(result, 3, 0);
      ELEM(Jr, jj, 1) = -2 * d2 * ELEM(result, 4, 0) * ELEM(result, 4, 0);
      ELEM(Jr, jj, 2) = -2 * d3 * ELEM(result, 5, 0) * ELEM(result, 5, 0);
      ELEM(Jr, jj, 3) = 2 * ELEM(result, 3, 0) * d1 * d1;
      ELEM(Jr, jj, 4) = 2 * ELEM(result, 4, 0) * d2 * d2;
      ELEM(Jr, jj, 5) = 2 * ELEM(result, 5, 0) * d3 * d3;
      t1 = ELEM(result, 3, 0) * d1;
      t2 = ELEM(result, 4, 0) * d2;
      t3 = ELEM(result, 5, 0) * d3;
      ELEM(res, jj, 0) = t1 * t1 + t2 * t2 + t3 * t3 - k2;
    }
    matPseudo_inv(Jr, tmp1);
    matMult(tmp1, res, delta);
    matSub(result, delta, result);
    if (matNorm(delta) < tol) {
      return 1;
    }
  }
  return 1;
}

//------------------Quadratic form (sort of)----------------------//
// returns matrix C=A*B*(~A)

void QuadProd(Matrix A, Matrix B, Matrix result)
{
  int16_t i, j, n, ii;
  float tmp;
  matZeros(result);
  for (n = 0; n < A.rows; n++) {
    for (i = 0; i < A.cols; i++) {
      tmp = 0.0;
      for (j = 0; j < A.cols; j++) {
        tmp += ELEM(A, n, j) * ELEM(B, i, j);
      }
      for (ii = 0; ii < A.rows; ii++) {
        ELEM(result, ii, n) += ELEM(A, ii, i) * tmp;
      }
    }
  }
  return;
}
