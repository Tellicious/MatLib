//
//  Matrix.c
//
//
//  Created by Andrea Vivani on 11/08/16.
//  Copyright © 2016 Andrea Vivani. All rights reserved.
//

#include "Matrix.h"
#include "math.h"
#include "NumMethods.h"
#include <stdlib.h>
#include <string.h>

//==========================================Assignment=============================================//

//-----------------------Constructor with external data-----------------------//
Matrix newMatrix(uint8_t rows, uint8_t cols, float *data)
{
  Matrix matrix;
  matrix.rows = rows;
  matrix.cols = cols;
  matrix.data = data;
  return matrix;
}

//---------------------Identity Matrix----------------------//
void matIdentity(Matrix matrix)
{
  uint16_t ii;
    matZeros(matrix);
    for (ii = 0; ii < ((matrix.cols < matrix.rows)? matrix.cols : matrix.rows); ii++)
        ELEM(matrix, ii, ii) = 1.0f;
  return;
}

//----------------------Zeros Matrix-----------------------//
void matZeros(Matrix matrix)
{
  uint16_t ii;
  for (ii = 0; ii < (matrix.cols * matrix.rows); ii++)
    matrix.data[ii] = 0.0f;
  return;
}
//==========================================Operations=============================================//
//--------------------Matrix copy data----------------------//
void matCopyData(Matrix input, Matrix output){
  uint16_t ii;
  for (ii = 0; ii < (input.cols * input.rows); ii++) {
    output.data[ii] = input.data[ii];
  }
  return;
}

//--------------------Matrix addition----------------------//
void matAdd(Matrix lhs, Matrix rhs, Matrix result)
{
  uint16_t ii;
  for (ii = 0; ii < (lhs.cols * lhs.rows); ii++) {
    result.data[ii] = lhs.data[ii] + rhs.data[ii];
  }
  return;
}

//--------------------Scalar addition----------------------//
void matAddScalar(Matrix lhs, float sc, Matrix result)
{
  uint16_t ii;
  for (ii = 0; ii < (lhs.cols * lhs.rows); ii++) {
    result.data[ii] = lhs.data[ii] + sc;
  }
  return;
}

//------------------Matrix subtraction--------------------//
void matSub(Matrix lhs, Matrix rhs, Matrix result)
{
  uint16_t ii;
  for (ii = 0; ii < (lhs.cols * lhs.rows); ii++) {
    result.data[ii] = lhs.data[ii] - rhs.data[ii];
  }
  return;
}

//---------------Matrix multiplication------------------//
uint8_t matMult(Matrix lhs, Matrix rhs, Matrix result)
{
  uint8_t i, j, k;
  if (lhs.cols != rhs.rows) {
    return 0;
  }
  matZeros(result);
  for (i = 0; i < lhs.rows; i++)
    for (j = 0; j < rhs.cols; j++)
      for (k = 0; k < lhs.cols; k++)
        ELEM(result, i, j) += ELEM(lhs, i, k) * ELEM(rhs, k, j);
  return 1;
}

//------Matrix multiplication with lhs transposed------//
uint8_t matMult_lhsT(Matrix lhs, Matrix rhs, Matrix result)
{
  uint8_t i, j, k;
  if (lhs.rows != rhs.rows) {
    return 0;
  }
  matZeros(result);
  for (i = 0; i < lhs.cols; i++)
    for (j = 0; j < rhs.cols; j++)
      for (k = 0; k < lhs.rows; k++)
        ELEM(result, i, j) += ELEM(lhs, k, i) * ELEM(rhs, k, j);
  return 1;
}

//------Matrix multiplication with rhs transposed------//
uint8_t matMult_rhsT(Matrix lhs, Matrix rhs, Matrix result)
{
  uint8_t i, j, k;
  if (lhs.cols != rhs.cols) {
    return 0;
  }
  matZeros(result);
  for (i = 0; i < lhs.rows; i++)
    for (j = 0; j < rhs.rows; j++)
      for (k = 0; k < lhs.cols; k++)
        ELEM(result, i, j) += ELEM(lhs, i, k) * ELEM(rhs, j, k);
  return 1;
}

//---------------Scalar multiplication------------------//
void matMultScalar(Matrix lhs, float sc, Matrix result)
{
  uint16_t ii;
  for (ii = 0; ii < (lhs.cols * lhs.rows); ii++) {
    result.data[ii] = lhs.data[ii] * sc;
  }
  return;
}

//--------------------Inverse LU------------------------//
void matInversed(Matrix lhs, Matrix result)
{
  float _eye_data[lhs.rows * lhs.cols];
  Matrix Eye = newMatrix(lhs.rows, lhs.cols, _eye_data);
  matIdentity(Eye);
  LinSolveLU(lhs, Eye, result);
  return;
}

//-----------------Robust Inverse LUP-------------------//
void matInversed_rob(Matrix lhs, Matrix result)
{
  float _eye_data[lhs.rows * lhs.cols];
  Matrix Eye = newMatrix(lhs.rows, lhs.cols, _eye_data);
  matIdentity(Eye);
  LinSolveLUP(lhs, Eye, result);
  return;
}

//-----------------Transposed--------------------//
void matTrans(Matrix lhs, Matrix result)
{
  uint8_t ii, jj;
  for (ii = 0; ii < lhs.rows; ii++)
    for (jj = 0; jj < lhs.cols; jj++)
      ELEM(result, jj, ii) = ELEM(lhs, ii, jj);
  return;
}

//-----------------Nomalized--------------------//
void matNormalized(Matrix lhs, Matrix result)
{
  float k = 1.0f / matNorm(lhs);
  matMultScalar(lhs, k, result);
  return;
}

//-------Moore-Penrose pseudo inverse---------//
void matPseudo_inv(Matrix lhs, Matrix result)
{
  float _tran_data[lhs.cols * lhs.rows];
  float _mult1_data[lhs.cols * lhs.cols];
  Matrix tran = newMatrix(lhs.cols, lhs.rows, _tran_data);
  Matrix mult1 = newMatrix(lhs.cols, lhs.cols, _mult1_data);
  matTrans(lhs, tran);
  matMult(tran, lhs, mult1);
  LinSolveLU(mult1, tran, result);
  return;
}

//=======================================Matrix Data=========================================//
//---------Returns one single element---------//
float matGet(Matrix matrix, uint8_t i, uint8_t j)
{
  return ELEM(matrix, i, j);
}

//-----------Returns the determinant----------//
float matDet(Matrix matrix)
{
  float _L_data[matrix.rows * matrix.rows];
  float _U_data[matrix.rows * matrix.rows];
  float _P_data[matrix.rows];
  Matrix L = newMatrix(matrix.rows, matrix.rows, _L_data);
  Matrix U = newMatrix(matrix.rows, matrix.rows, _U_data);
  Matrix P = newMatrix(matrix.rows, 1, _P_data);
  int16_t ii;
  int8_t det_f;
  float determinant = 1.0f;

  if (matrix.rows != matrix.cols) {
    return 0.0f;
  }

  if (LU_Cormen(matrix, L, U)) {
    for (ii = 0; ii < matrix.rows; ii++) {
      determinant *= ELEM(U, ii, ii);
    }
  }

  else {
    det_f = LUP_Cormen(matrix, L, U, P);
    if (det_f) {
      for (ii = 0; ii < matrix.rows; ii++) {
        determinant *= ELEM(U, ii, ii);
      }
      determinant *= det_f;
    }
    else {
      determinant = 0.0f;
    }
  }

  return determinant;
}

//-------------Returns the norm--------------//
float matNorm(Matrix matrix)
{
  float result = 0.0f;
  uint16_t i;
  for (i = 0; i < (matrix.rows * matrix.cols); i++) {
    result += matrix.data[i] * matrix.data[i];
  }
  result = sqrtf(result);
  return result;
}
