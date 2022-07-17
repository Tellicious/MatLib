//
//  Matrix.h
//
//
//  Created by Andrea Vivani on 11/08/16.
//  Copyright © 2016 Andrea Vivani. All rights reserved.
//

#ifndef Matrix_h
#define Matrix_h
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define ELEM(m, i, j) (m->data[(i) * m->cols + (j)])

struct _Matrix {
  uint8_t rows;
  uint8_t cols;
  float *data;
};
typedef struct _Matrix *Matrix;

//==================Assignment===================//
Matrix newMatrix(uint8_t rows, uint8_t cols);
Matrix newMatrixData(uint8_t rows, uint8_t cols, float *data); //creates a new matrix with the provided data
Matrix copyMatrix(Matrix matrix); //creates a new matrix copying the provided one
void matResize(Matrix matrix, uint8_t rows, uint8_t cols); //resize existing matrix
void matFree(Matrix matrix); //destructor
void matIdentity(Matrix matrix); // sets the matrix as an identity matrix
void matZeros(Matrix matrix); // zeros the matrix
//==================Operations==================//
void matCopyData(Matrix input, Matrix output); //copy data of Matrix input to Matrix output
void matAdd(Matrix lhs, Matrix rhs, Matrix result); // Matrix addition
void matAddScalar(Matrix lhs, float sc, Matrix result); // Matrix and scalar addition
void matSub(Matrix lhs, Matrix rhs, Matrix result); // Matrix subtraction
uint8_t matMult(Matrix lhs, Matrix rhs, Matrix result); // Matrix multiplication
void matMultScalar(Matrix lhs, float sc, Matrix result); // Matrix scalar multiplication
void matInversed(Matrix lhs, Matrix result); // inversed Matrix
void matInversed_rob(Matrix lhs, Matrix result); //robust inversed Matrix using LUP decomposition
void matTrans(Matrix lhs, Matrix result); //transposed Matrix
void matNormalized(Matrix lhs, Matrix result); //normalized Matrix
void matPseudo_inv(Matrix lhs, Matrix result); //Moore-Penrose pseudo inverse
//==================Matrix Data==================//
float matGet(Matrix matrix, uint8_t i, uint8_t j); //retrieves one single element
float matDet(Matrix matrix); //retrieves the determinant of the Matrix
float matNorm(Matrix matrix); //retrives the norm of the Matrix

#ifdef __cplusplus
}
#endif

#endif /* Matrix_h */
