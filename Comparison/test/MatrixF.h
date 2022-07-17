//
//  MatrixF.h
//
//
//  Created by Andrea Vivani on 11/08/16.
//  Copyright © 2016 Andrea Vivani. All rights reserved.
//

#ifndef MatrixF_h
#define MatrixF_h
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define MAT(m, i, j) (m.mat[(i) * m.cols + (j)])
    
typedef struct _MatrixF{
    uint8_t rows;
    uint8_t cols;
    float mat[200];
} MatrixF;
//typedef struct _MatrixF* MatrixF;

//==================Assignment===================//
MatrixF newMatrix(uint8_t rows, uint8_t cols);
MatrixF newMatrixData(uint8_t rows, uint8_t cols, float* data); //creates a new matrix with the provided data
MatrixF matIdentity(uint8_t rows, uint8_t cols); // returns an identity matrix
MatrixF matZeros(uint8_t rows, uint8_t cols); // returns a zeroed matrix
//==================Operations==================//
MatrixF matAdd(MatrixF lhs, MatrixF rhs); // Matrix addition
MatrixF matAddScalar(MatrixF lhs, float sc); // Matrix and scalar addition
MatrixF matSub(MatrixF lhs, MatrixF rhs); // Matrix subtraction
MatrixF matSubScalar(MatrixF lhs, float sc); // Matrix and scalar subtraction
MatrixF matMult(MatrixF lhs, MatrixF rhs);// Matrix multiplication
MatrixF matMultScalar(MatrixF lhs, float sc);// Matrix scalar multiplication
MatrixF matInversed(MatrixF lhs); // inversed Matrix
MatrixF matInversed_rob(MatrixF lhs); //robust inversed Matrix using LUP decomposition
MatrixF matTrans(MatrixF lhs); //transposed Matrix
MatrixF matNormalized(MatrixF lhs); //normalized Matrix
MatrixF matPseudo_inv(MatrixF lhs); //Moore-Penrose pseudo inverse
//==================Matrix Data==================//
float matGet(MatrixF matrix, uint8_t i, uint8_t j); //retrieves one single element
float matDet(MatrixF matrix); //retrieves the determinant of the Matrix
float matNorm(MatrixF matrix); //retrives the norm of the Matrix

#ifdef __cplusplus
}
#endif

#endif /* MatrixF_h */
