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
//-----------------------Constructor------------------------//
Matrix newMatrix(uint8_t rows, uint8_t cols)
{
  Matrix matrix = malloc(sizeof(struct _Matrix));
  matrix->rows = rows;
  matrix->cols = cols;
  matrix->data = calloc(rows * cols, sizeof(float));
  return matrix;
}

//-----------------------Constructor with provided data-----------------------//
Matrix newMatrixData(uint8_t rows, uint8_t cols, float *data)
{
  Matrix matrix = newMatrix(rows, cols);
  uint16_t ii;
  for (ii = 0; ii < (cols * rows); ii++)
    matrix->data[ii] = data[ii];
  return matrix;
}

//----------------------------Copy other matrix-------------------------------//
Matrix copyMatrix(Matrix matrix)
{
  Matrix copy = NULL;
  if (!matrix)
    return NULL;
  // create a new matrix to hold the copy
  copy = newMatrix(matrix->rows, matrix->cols);
  if (matrix->data)
    memcpy(copy->data, matrix->data, matrix->rows * matrix->cols * sizeof(float));
  return copy;
}

//-----------------------------Resize matrix----------------------------------//
void matResize(Matrix matrix, uint8_t rows, uint8_t cols)
{
  matrix->rows = rows;
  matrix->cols = cols;
  matrix->data = realloc(matrix->data, rows * cols * sizeof(float));
  return;
}

//---------------------------Destructor--------------------------//
void matFree(Matrix matrix)
{
  if (matrix->data)
    free(matrix->data);
  free(matrix);
}

//---------------------Identity Matrix----------------------//
void matIdentity(Matrix matrix)
{
  uint16_t ii;
  for (ii = 0; ii < (matrix->cols * matrix->rows); ii++)
    matrix->data[ii] = (
        ((ii / matrix->cols) == (ii % matrix->cols)) ? 1.0f : 0.0f);
  return;
}

//----------------------Zeros Matrix-----------------------//
void matZeros(Matrix matrix)
{
  uint16_t ii;
  for (ii = 0; ii < (matrix->cols * matrix->rows); ii++)
    matrix->data[ii] = 0.0f;
  return;
}
//==========================================Operations=============================================//
//--------------------Matrix copy data----------------------//
void matCopyData(Matrix input, Matrix output){
  uint16_t ii;
  for (ii = 0; ii < (input->cols * input->rows); ii++) {
    output->data[ii] = input->data[ii];
  }
  return;
}

//--------------------Matrix addition----------------------//
void matAdd(Matrix lhs, Matrix rhs, Matrix result)
{
  uint16_t ii;
  matResize(result, lhs->rows, lhs->cols);
  for (ii = 0; ii < (lhs->cols * lhs->rows); ii++) {
    result->data[ii] = lhs->data[ii] + rhs->data[ii];
  }
  return;
}

//--------------------Scalar addition----------------------//
void matAddScalar(Matrix lhs, float sc, Matrix result)
{
  uint16_t ii;
  matResize(result, lhs->rows, lhs->cols);
  for (ii = 0; ii < (lhs->cols * lhs->rows); ii++) {
    result->data[ii] = lhs->data[ii] + sc;
  }
  return;
}

//------------------Matrix subtraction--------------------//
void matSub(Matrix lhs, Matrix rhs, Matrix result)
{
  uint16_t ii;
  matResize(result, lhs->rows, lhs->cols);
  for (ii = 0; ii < (lhs->cols * lhs->rows); ii++) {
    result->data[ii] = lhs->data[ii] - rhs->data[ii];
  }
  return;
}

//---------------Matrix multiplication------------------//
uint8_t matMult(Matrix lhs, Matrix rhs, Matrix result)
{
  uint8_t i, j, k;
  matResize(result, lhs->rows, rhs->cols);
  if (lhs->cols != rhs->rows) {
    return 0;
  }
  matZeros(result);
  for (i = 0; i < lhs->rows; i++)
    for (j = 0; j < rhs->cols; j++)
      for (k = 0; k < lhs->cols; k++)
        ELEM(result, i, j) += ELEM(lhs, i, k) * ELEM(rhs, k, j);
  return 1;
}

//---------------Scalar multiplication------------------//
void matMultScalar(Matrix lhs, float sc, Matrix result)
{
  uint16_t ii;
  matResize(result, lhs->rows, lhs->cols);
  for (ii = 0; ii < (lhs->cols * lhs->rows); ii++) {
    result->data[ii] = lhs->data[ii] * sc;
  }
  return;
}

//--------------------Inverse LU------------------------//
void matInversed(Matrix lhs, Matrix result)
{
  Matrix Eye = newMatrix(lhs->rows, lhs->cols);
  matResize(result, lhs->rows, lhs->cols);
  matIdentity(Eye);
  LinSolveLU(lhs, Eye, result);
  matFree(Eye);
  return;
}

//-----------------Robust Inverse LUP-------------------//
void matInversed_rob(Matrix lhs, Matrix result)
{
  Matrix Eye = newMatrix(lhs->rows, lhs->cols);
  matResize(result, lhs->rows, lhs->cols);
  matIdentity(Eye);
  LinSolveLUP(lhs, Eye, result);
  matFree(Eye);
  return;
}

//-----------------Transposed--------------------//
void matTrans(Matrix lhs, Matrix result)
{
  uint8_t ii, jj;
  matResize(result, lhs->cols, lhs->rows);
  for (ii = 0; ii < lhs->rows; ii++)
    for (jj = 0; jj < lhs->cols; jj++)
      ELEM(result, jj, ii) = ELEM(lhs, ii, jj);
  return;
}

//-----------------Nomalized--------------------//
void matNormalized(Matrix lhs, Matrix result)
{
  float k = 1.0f / matNorm(lhs);
  matResize(result, lhs->rows, lhs->cols);
  matMultScalar(lhs, k, result);
  return;
}

//-------Moore-Penrose pseudo inverse---------//
void matPseudo_inv(Matrix lhs, Matrix result)
{
  Matrix tran = newMatrix(lhs->cols, lhs->rows);
  Matrix mult1 = newMatrix(lhs->cols, lhs->cols);
  matResize(result, lhs->cols, lhs->rows);
  matTrans(lhs, tran);
  matMult(tran, lhs, mult1);
  LinSolveLU(mult1, tran, result);
  matFree(tran);
  matFree(mult1);
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
  Matrix L = newMatrix(matrix->rows, matrix->rows);
  Matrix U = newMatrix(matrix->rows, matrix->rows);
  Matrix P = newMatrix(matrix->rows, 1);
  int16_t ii;
  int8_t det_f;
  float determinant = 1.0f;

  if (matrix->rows != matrix->cols) {
    matFree(L);
    matFree(U);
    matFree(P);
    return 0.0f;
  }

  if (LU_Cormen(matrix, L, U)) {
    for (ii = 0; ii < matrix->rows; ii++) {
      determinant *= ELEM(U, ii, ii);
    }
  }

  else {
    det_f = LUP_Cormen(matrix, L, U, P);
    if (det_f) {
      for (ii = 0; ii < matrix->rows; ii++) {
        determinant *= ELEM(U, ii, ii);
      }
      determinant *= det_f;
    }
    else {
      determinant = 0.0f;
    }
  }

  matFree(L);
  matFree(U);
  matFree(P);
  return determinant;
}

//-------------Returns the norm--------------//
float matNorm(Matrix matrix)
{
  float result = 0.0f;
  uint16_t i;
  for (i = 0; i < (matrix->rows * matrix->cols); i++) {
    result += matrix->data[i] * matrix->data[i];
  }
  result = sqrtf(result);
  return result;
}
