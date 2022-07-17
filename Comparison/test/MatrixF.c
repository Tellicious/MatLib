//
//  MatrixF.c
//
//
//  Created by Andrea Vivani on 11/08/16.
//  Copyright © 2016 Andrea Vivani. All rights reserved.
//

#include "MatrixF.h"
#include "math.h"
#include "NumMethodsF.h"
#include <stdlib.h>

//==========================================Assignment=============================================//
//-----------------------Constructor------------------------//
MatrixF newMatrix(uint8_t rows, uint8_t cols){
    MatrixF matrix; //= (MatrixF) malloc(sizeof(struct _MatrixF));
    matrix.rows = rows;
    matrix.cols = cols;
    //matrix.mat = (float*) malloc(rows * cols * sizeof(float));
    return matrix;
}

//-----------------------Constructor with provided data-----------------------//
MatrixF newMatrixData(uint8_t rows, uint8_t cols, float* data){
    MatrixF matrix = newMatrix(rows, cols);
    uint16_t ii;
    for (ii = 0; ii < (cols * rows); ii++)
        matrix.mat[ii] = data[ii];
    return matrix;
}

//---------------------Identity Matrix----------------------//
MatrixF matIdentity(uint8_t rows, uint8_t cols){
    MatrixF matrix = newMatrix(rows, cols);
    uint16_t ii;
    for (ii = 0; ii < (cols * rows); ii++)
        matrix.mat[ii] = (((ii / cols) == (ii % cols))? 1.0f : 0.0f);
    return matrix;
}

//----------------------Zeros Matrix-----------------------//
MatrixF matZeros(uint8_t rows, uint8_t cols){
    MatrixF matrix = newMatrix(rows, cols);
    uint16_t ii;
    for (ii = 0; ii < (cols * rows); ii++)
        matrix.mat[ii] = 0.0f;
    return matrix;
}

//==========================================Operations=============================================//
//--------------------Matrix addition----------------------//
MatrixF matAdd(MatrixF lhs, MatrixF rhs){
    MatrixF m = newMatrix(lhs.rows, lhs.cols);
    uint16_t ii;
    for (ii = 0; ii < (lhs.cols * lhs.rows); ii++) {
        m.mat[ii] = lhs.mat[ii] + rhs.mat[ii];
    }
    return m;
}

//--------------------Scalar addition----------------------//
MatrixF matAddScalar(MatrixF lhs, float sc){
    MatrixF m = newMatrix(lhs.rows, lhs.cols);
    uint8_t ii;
    for (ii = 0; ii < (lhs.cols * lhs.rows); ii++) {
        m.mat[ii] = lhs.mat[ii] + sc;
    }
    return m;
}

//------------------Matrix subtraction--------------------//
MatrixF matSub(MatrixF lhs, MatrixF rhs){
    MatrixF m = newMatrix(lhs.rows, lhs.cols);
    uint16_t ii;
    for (ii = 0; ii < (lhs.cols * lhs.rows); ii++) {
        m.mat[ii] = lhs.mat[ii] - rhs.mat[ii];
    }
    return m;
}

//------------------Scalar subtraction--------------------//
MatrixF matSubScalar(MatrixF lhs, float sc){
    MatrixF m = newMatrix(lhs.rows, lhs.cols);
    uint16_t ii;
    for (ii = 0; ii < (lhs.cols * lhs.rows); ii++) {
        m.mat[ii] = lhs.mat[ii] - sc;
    }
    return m;
}

//---------------Matrix multiplication------------------//
MatrixF matMult(MatrixF lhs, MatrixF rhs){
    uint8_t i, j, k;
    if (lhs.cols != rhs.rows) {
        return matZeros(1, 1);
    }
    MatrixF result = matZeros(lhs.rows, rhs.cols);
    for (i=0; i < lhs.rows; i++)
        for (j=0; j < rhs.cols; j++)
            for (k=0; k < lhs.cols; k++)
                MAT(result, i, j) += MAT(lhs, i, k) * MAT(rhs, k, j);
    return result;
}

//---------------Scalar multiplication------------------//
MatrixF matMultScalar(MatrixF lhs, float sc){
    MatrixF m = newMatrix(lhs.rows, lhs.cols);
    uint16_t ii;
    for (ii = 0; ii < (lhs.cols * lhs.rows); ii++) {
        m.mat[ii] = lhs.mat[ii] * sc;
    }
    return m;
}

//--------------------Inverse LU------------------------//
MatrixF matInversed(MatrixF lhs);
//-----------------Robust Inverse LUP-------------------//
MatrixF matInversed_rob(MatrixF lhs);
//-----------------Transposed--------------------//
MatrixF matTrans(MatrixF lhs){
    MatrixF m = newMatrix(lhs.cols, lhs.rows);
    uint8_t ii, jj;
    for (ii = 0; ii < lhs.rows; ii++)
        for (jj = 0; jj < lhs.cols; jj++)
            MAT(m, jj, ii) = MAT(lhs, ii, jj);
    return m;
}

//-----------------Nomalized--------------------//
MatrixF matNormalized(MatrixF lhs){
    float k = 1.0f / matNorm(lhs);
    return matMultScalar(lhs, k);
}

//-------Moore-Penrose pseudo inverse---------//
MatrixF matPseudo_inv(MatrixF lhs) {
    MatrixF tran = matTrans(lhs);
    return LinSolveLUF(matMult(tran, lhs), tran);
}

//=======================================Matrix Data=========================================//
//---------Returns one single element---------//
float matGet(MatrixF matrix, uint8_t i, uint8_t j){
    return MAT(matrix, i, j);
}

//-----------Returns the determinant----------//
float matDet(MatrixF matrix);

//-------------Returns the norm--------------//
float matNorm(MatrixF matrix){
    float result = 0.0f;
    uint16_t i;
    for (i = 0; i < (matrix.rows * matrix.cols); i++){
        result += matrix.mat[i] * matrix.mat[i];
    }
    result = sqrtf(result);
    return result;
}
