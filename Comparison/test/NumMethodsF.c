//
//  NumMethodsF.c
//
//
//  Created by Andrea Vivani on 11/08/16.
//  Copyright © 2016 Andrea Vivani. All rights reserved.
//

#include "NumMethodsF.h"
#include "math.h"
#define MATP(m, i, j) (m->mat[(i) * m->cols + (j)])

//-------------------Forward substitution----------------------//
// assumes that the matrix A is already a lower triangular one. No check!

MatrixF fwsubF(MatrixF A, MatrixF B){
    MatrixF result = newMatrix(A.cols, B.cols);
    int16_t i, j, k;
    float tmp;
    for (k = 0; k < B.cols; k++){
        MAT(result, 0, k) = MAT(B, 0, k) / MAT(A, 0, 0);
        for (i = 1; i < A.rows; i++){
            tmp = 0.0;
            for (j = 0; j < i; j++){
                tmp += MAT(A, i, j) * (MAT(result, j, k));
            }
            MAT(result, i, k) = (MAT(B, i, k) - tmp) / MAT(A, i, i);
        }
    }
    return result;
};

//---------------Forward substitution with permutation-------------------//
// assumes that the matrix A is already a lower triangular one. No check!

MatrixF fwsubPermF(MatrixF A, MatrixF B, MatrixF P){
    MatrixF result = newMatrix(A.cols, B.cols);
    int16_t i, j, k;
    float tmp;
    for (k = 0; k < B.cols; k++){
        MAT(result, 0, k) = MAT(B, (uint8_t) MAT(P, 0, 0), k) / MAT(A, 0, 0);
        for (i = 1; i < A.rows; i++){
            tmp = 0.0;
            for (j = 0; j < i; j++){
                tmp += MAT(A, i, j) * MAT(result, j, k);
            }
            MAT(result, i, k)= (MAT(B, (uint8_t) MAT(P, i, 0), k) - tmp) / MAT(A, i, i);
        }
    }
    return result;
};

//-------------------Backward substitution----------------------//
// assumes that the matrix A is already an upper triangular one. No check!

MatrixF bksubF(MatrixF A, MatrixF B){
    uint8_t ncolsA = A.cols;
    MatrixF result = newMatrix(ncolsA, B.cols);
    int16_t i, j, k;
    float tmp;
    for (k = 0; k < B.cols; k++){
        MAT(result, ncolsA - 1, k) = MAT(B, ncolsA - 1, k) / MAT(A, ncolsA - 1, ncolsA - 1);
        for (i = A.rows - 2; i >= 0; i--){
            tmp = 0.0;
            for (j = ncolsA - 1; j > i; j--){
                tmp += MAT(A, i, j) * MAT(result, j, k);
            }
            MAT(result, i, k)= (MAT(B, i, k) - tmp) / MAT(A, i, i);
        }
    }
    return result;
};

//--------------Backward substitution with permutation-----------------//
// assumes that the matrix A is already an upper triangular one. No check!

MatrixF bksubPermF(MatrixF A, MatrixF B, MatrixF P){
    uint8_t ncolsA = A.cols;
    MatrixF result = newMatrix(ncolsA, B.cols);
    int16_t i, j, k;
    float tmp;
    for (k = 0; k < B.cols; k++){
        MAT(result, ncolsA - 1, k) = MAT(B, (uint8_t) MAT(P, ncolsA - 1, 0), k) / MAT(A, ncolsA - 1,ncolsA - 1);
        for (i = A.rows - 2; i >= 0; i--){
            tmp=0.0;
            for (j = ncolsA - 1; j > i; j--){
                tmp += MAT(A, i, j) * MAT(result, j, k);
            }
            MAT(result, i, k)= (MAT(B, (uint8_t) MAT(P, i, 0), k) - tmp) / MAT(A, i, i);
        }
    }
    return result;
};

//-------------------------LU factorization using Crout's Method--------------------------------//
// factorizes the A matrix as the product of a unit upper triangular matrix U and a lower triangular matrix L

uint8_t LU_CroutF(MatrixF A, MatrixF* L, MatrixF* U){
    int16_t ii, jj, kk;
    float sum = 0.0;
    uint8_t nrowsA = A.rows;
    *U = matIdentity(U->rows, U->cols);
    *L = matZeros(L->rows, L->cols);
    for (jj = 0; jj < nrowsA; jj++) {
        for (ii = jj; ii < nrowsA; ii++) {
            sum = 0.0;
            for (kk = 0; kk < jj; kk++) {
                sum += MATP(L, ii, kk) * MATP(U, kk, jj);
            }
            MATP(L, ii, jj) = MAT(A, ii, jj) - sum;
        }
        
        for (ii = jj; ii < nrowsA; ii++) {
            sum = 0;
            for(kk = 0; kk < jj; kk++) {
                sum += MATP(L, jj, kk) * MATP(U, kk, ii);
            }
            if (MATP(L, jj, jj) == 0) {
                return 0;
            }
            MATP(U, jj, ii) = (MAT(A, jj, ii) - sum) / MATP(L, jj, jj);
        }
    }
    return 1;
};

//-------------------------LU factorization using Cormen's Method--------------------------------//
// factorizes the A matrix as the product of a unit upper triangular matrix U and a lower triangular matrix L

uint8_t LU_CormenF(MatrixF A, MatrixF* L, MatrixF* U){
    int16_t i, j, k;
    float tmp;
    uint8_t nrowsA = A.rows;
    *U = matZeros(U->rows, U->cols);
    *L = matIdentity(L->rows, L->cols);
    
    for (k = 0; k < nrowsA; k++) {
        MATP(U, k, k) = MAT(A, k, k);
        if (MAT(A, k,k) == 0){
            return 0;
        }
        tmp = 1.0 / MATP(U, k, k);
        for (i = k + 1; i < nrowsA; i++) {
            MATP(L, i, k) = MAT(A, i, k) * tmp;
            MATP(U, k, i) = MAT(A, k, i);
        }
        for (i = k + 1; i < nrowsA; i++) {
            for (j = k + 1; j < nrowsA; j++) {
                MAT(A, i, j) -= MATP(L, i, k) * MATP(U, k, j);
            }
        }
    }
    return 1;
};

//-----------------------LUP factorization using Cormen's Method------------------------------//
// factorizes the A matrix as the product of a upper triangular matrix U and a unit lower triangular matrix L
// returns the factor that has to be multiplied to the determinant of U in order to obtain the correct value

int8_t LUP_CormenF(MatrixF A, MatrixF* L, MatrixF* U, MatrixF* P){
    int16_t i, j, k;
    float tmp, tmp2;
    uint8_t nrowsA = A.rows;
    int16_t pivrow;
    int8_t d_mult=1; // determinant multiplying factor
    *U = matZeros(U->rows, U->cols);
    *L = matIdentity(L->rows, L->cols);
    // initialization
    for (i = 0; i < nrowsA; i++){
        MATP(P, i, 0) = i;
    }
    
    // outer loop over diagonal pivots
    for (k = 0; k < nrowsA-1; k++) {
        
        // inner loop to find the largest pivot
        pivrow = k;
        tmp = fabs(MAT(A, k, k));
        for (i = k + 1 ; i < nrowsA; i++){
            tmp2 = fabs(MAT(A, i, k));
            if (tmp2 > tmp){
                tmp=tmp2;
                pivrow = i;
            }
        }
        // check for singularity
        if (MAT(A, pivrow, k)==0) {
            return 0;
        }
        
        // swap rows
        if (pivrow != k) {
            tmp = MATP(P, k, 0);
            MATP(P, k, 0) = MATP(P, pivrow, 0);
            MATP(P, pivrow, 0) = tmp;
            d_mult *= -1;
            
            for (j = 0; j < nrowsA; j++){
                tmp = MAT(A, k, j);
                MAT(A, k, j) = MAT(A, pivrow, j);
                MAT(A, pivrow, j) = tmp;
            }
        }
        tmp = 1.0 / MAT(A, k, k);
        // Gaussian elimination
        for (i = k + 1; i < nrowsA; i++) { // iterate down rows
            MAT(A, i, k) *= tmp;
            for (j = k + 1; j < nrowsA; j++){ // iterate across rows
                MAT(A, i, j) -= MAT(A, i, k) * MAT(A, k, j);
            }
        }
    }
    for (k = 0; k < nrowsA; k++){
        MATP(U, k, k) = MAT(A, k, k);
        for (j = k + 1; j < nrowsA; j++){
            MATP(L, j, k) = MAT(A, j, k);
            MATP(U, k, j) = MAT(A, k, j);
        }
    }
    return d_mult;
};

//-----------------------Linear system solver using LU factorization---------------------------//
// solves the linear system A*X=B, where A is a n-by-n matrix and B an n-by-m matrix, giving the n-by-m matrix X

MatrixF LinSolveLUF(MatrixF A, MatrixF B) {
    MatrixF L = matZeros(A.rows, A.cols);
    MatrixF U = matZeros(A.rows, A.cols);
    LU_CormenF(A, &L, &U);
    return bksubF(U, fwsubF(L, B));
};

//----------------------Linear system solver using LUP factorization--------------------------//
// solves the linear system A*X=B, where A is a n-by-n matrix and B an n-by-m matrix, giving the n-by-m matrix X

MatrixF LinSolveLUPF(MatrixF A, MatrixF B) {
    MatrixF L = matZeros(A.rows, A.cols);
    MatrixF U = matZeros(A.rows, A.cols);
    MatrixF P = matZeros(A.rows, 1);
    LUP_CormenF(A, &L, &U, &P);
    return bksubF(U, fwsubPermF(L, B, P));
};

//------------Linear system solver using Gauss elimination with partial pivoting---------------//
// solves the linear system A*X=B, where A is a n-by-n matrix and B an n-by-m matrix, giving the n-by-m matrix X

MatrixF LinSolveGaussF(MatrixF A, MatrixF B) {
    uint8_t pivrow = 0;     // keeps track of current pivot row
    uint8_t ncolsA = A.cols;
    uint8_t ncolsB = B.cols;
    uint8_t k, i, j;      // k: overall index along diagonals; i: row index; j: col index
    float tmp;      // used for finding max value and making row swaps
    float tmp2;  // used to store abs when finding max value and to store coefficient value when eliminating values
    
    for (k = 0; k < (ncolsA-1); k++){
        
        // find pivot row, the row with biggest entry in current column
        tmp = fabs(MAT(A, k, k));
        pivrow = k;
        for (i = k + 1; i < ncolsA; i++){
            tmp2 = fabs(MAT(A, i, k)); // 'Avoid using other functions inside abs()?'
            if (tmp2 > tmp){
                tmp = tmp2;
                pivrow = i;
            }
        }
        
        // check for singular Matrix
        if (MAT(A, pivrow, k) == 0.0){
            return matZeros(1, 1);
        }
        
        // Execute pivot (row swap) if needed
        if (pivrow != k){
            // swap row k of matrix A with pivrow
            for (j = k; j < ncolsA; j++){
                tmp = MAT(A, k, j);
                MAT(A, k, j) = MAT(A, pivrow, j);
                MAT(A, pivrow, j) = tmp;
            }
            // swap row k of matrix B with pivrow
            for (j = 0; j < ncolsB; j++)
            {
                tmp = MAT(B, k, j);
                MAT(B, k, j) = MAT(B, pivrow, j);
                MAT(B, pivrow, j) = tmp;
            }
        }
        
        // Row reduction
        tmp = 1.0 / MAT(A, k, k);    // invert pivot element
        for (i = k + 1; i < ncolsA; i++){ // along rows
            tmp2 = MAT(A, i, k) * tmp;
            // Perform row reduction of A
            for (j = k + 1; j < ncolsA; j++){ //along columns of A
                MAT(A, i, j) -= tmp2 * MAT(A, k, j);
            }
            // Perform row reduction of B
            for (j = 0; j < ncolsB; j++){ //along columns of B
                MAT(B, i, j) -= tmp2 * MAT(B, k, j);
            }
        }
        
    }
    return(bksubF(A, B));
};

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

MatrixF GaussNewton_Sens_Cal_9F(MatrixF Data, float k, MatrixF X0, uint16_t nmax, float tol){
    MatrixF result = X0;
    uint8_t nrows = Data.rows;
    uint8_t ncols = Data.cols;
    float d1, d2, d3, rx1, rx2, rx3, t1, t2, t3;
    float k2 = k * k;
    uint16_t n_iter;
    uint8_t jj;
    MatrixF Jr = matZeros(nrows, 9);
    MatrixF res = matZeros(nrows, 1);
    MatrixF delta = newMatrix(9, 1);
    
    if ((nrows < 9) || (ncols != 3))
        return matZeros(1, 1);
    
    for (n_iter = 0; n_iter < nmax; n_iter++){
        for (jj = 0; jj < nrows; jj++){
            d1 = MAT(Data, jj, 0) - MAT(result, 0, 0);
            d2 = MAT(Data, jj, 1) - MAT(result, 1, 0);
            d3 = MAT(Data, jj, 2) - MAT(result, 2, 0);
            rx1 = -2 * (MAT(result, 3, 0) * d1 + MAT(result, 4, 0) * d2 + MAT(result, 5, 0) * d3);
            rx2 = -2 * (MAT(result, 4, 0) * d1 + MAT(result, 6, 0) * d2 + MAT(result, 7, 0) * d3);
            rx3 = -2 * (MAT(result, 5, 0) * d1 + MAT(result, 7, 0) * d2 + MAT(result, 8, 0) * d3);
            MAT(Jr, jj, 0) = MAT(result, 3, 0) * rx1 + MAT(result, 4, 0) * rx2 + MAT(result, 5, 0) * rx3;
            MAT(Jr, jj, 1) = MAT(result, 4, 0) * rx1 + MAT(result, 6, 0) * rx2 + MAT(result, 7, 0) * rx3;
            MAT(Jr, jj, 2) = MAT(result, 5, 0) * rx1 + MAT(result, 7, 0) * rx2 + MAT(result, 8, 0) * rx3;
            MAT(Jr, jj, 3) = -d1 * rx1;
            MAT(Jr, jj, 4) = -d2 * rx1 - d1 * rx2;
            MAT(Jr, jj, 5) = -d3 * rx1 - d1 * rx3;
            MAT(Jr, jj, 6) = -d2 * rx2;
            MAT(Jr, jj, 7) = -d3 * rx2 - d2 * rx3;
            MAT(Jr, jj, 8) = -d3 * rx3;
            t1 = MAT(result, 3, 0) * d1 + MAT(result, 4,0) * d2 + MAT(result, 5, 0) * d3;
            t2 = MAT(result, 4, 0) * d1 + MAT(result, 6,0) * d2 + MAT(result, 7, 0) * d3;
            t3 = MAT(result, 5, 0) * d1 + MAT(result, 7,0) * d2 + MAT(result, 8, 0) * d3;
            MAT(res, jj, 0) = t1 * t1 + t2 * t2 + t3 * t3 - k2;
        }

        delta = matMult(matPseudo_inv(Jr), res);
        result = matSub(result, delta);
        if (matNorm(delta) < tol)
            return result;
    }
    return result;
};


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

MatrixF GaussNewton_Sens_Cal_6F(MatrixF Data, float k, MatrixF X0, uint16_t nmax, float tol){
    MatrixF result = X0;
    uint8_t nrows = Data.rows;
    uint8_t ncols = Data.cols;
    float d1, d2, d3, t1, t2, t3;
    float k2 = k * k;
    uint16_t n_iter, jj;
    MatrixF Jr = matZeros(nrows, 6);
    MatrixF res = matZeros(nrows, 1);
    MatrixF delta = matZeros(6, 1);
    
    if ((nrows < 6) || (ncols != 3))
        return matZeros(1, 1);
    
    for (n_iter = 0; n_iter < nmax; n_iter++){
        for (jj = 0; jj < nrows; jj++){
            d1 = MAT(Data, jj, 0) - MAT(result, 0, 0);
            d2 = MAT(Data, jj, 1) - MAT(result, 1, 0);
            d3 = MAT(Data, jj, 2) - MAT(result, 2, 0);
            MAT(Jr, jj, 0) = -2 * d1 * MAT(result, 3, 0) * MAT(result, 3, 0);
            MAT(Jr, jj, 1) = -2 * d2 * MAT(result, 4, 0) * MAT(result, 4, 0);
            MAT(Jr, jj, 2) = -2 * d3 * MAT(result, 5, 0) * MAT(result, 5, 0);
            MAT(Jr, jj, 3) = 2 * MAT(result, 3, 0) * d1 * d1;
            MAT(Jr, jj, 4) = 2 * MAT(result, 4, 0) * d2 * d2;
            MAT(Jr, jj, 5) = 2 * MAT(result, 5, 0) * d3 * d3;
            t1 = MAT(result, 3, 0) * d1;
            t2 = MAT(result, 4, 0) * d2;
            t3 = MAT(result, 5, 0) * d3;
            MAT(res, jj, 0) = t1 * t1 + t2 * t2 + t3 * t3 - k2;
        }
        // pseudo-inverse*res (inv(Jr'*Jr)*Jr'*res)
        delta = matMult(matPseudo_inv(Jr), res);
        result = matSub(result, delta);
        if (matNorm(delta) < tol)
            return result;
    }
    return result;
};

//------------------Quadratic form (sort of)----------------------//
// returns matrix C=A*B*(~A)

MatrixF QuadProdF(MatrixF A, MatrixF B){
    MatrixF result = newMatrix(A.rows, A.rows);
    int16_t i, j, n, ii;
    float tmp;
    for (n = 0; n < A.rows; n++) {
        for (i = 0; i < A.cols; i++) {
            tmp = 0.0;
            for (j = 0; j < A.cols; j++){
                tmp += MAT(A, n, j) * MAT(B, i, j);
            }
            for (ii = 0; ii < A.rows; ii++) {
                MAT(result, ii, n) += MAT(A, ii, i) * tmp;
            }
        }
    }
    return result;
};
