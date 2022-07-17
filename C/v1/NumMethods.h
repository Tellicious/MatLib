//
//  NumMethods.h
//
//
//  Created by Andrea Vivani on 11/08/16.
//  Copyright © 2016 Andrea Vivani. All rights reserved.
//

#ifndef NumMethodsF_h
#define NumMethodsF_h
#include "Matrix.h"

#ifdef __cplusplus
extern "C" {
#endif
    //-------------------Forward substitution----------------------//
    // assumes that the matrix A is already a lower triangular one. No check!
    void fwsub(Matrix A, Matrix B, Matrix result);
    //---------------Forward substitution with permutation-------------------//
    // assumes that the matrix A is already a lower triangular one. No check!
    void fwsubPerm(Matrix A, Matrix B, Matrix P, Matrix result);
    //-------------------Backward substitution----------------------//
    // assumes that the matrix A is already an upper triangular one. No check!
    void bksub(Matrix A, Matrix B, Matrix result);
    //--------------Backward substitution with permutation-----------------//
    // assumes that the matrix A is already an upper triangular one. No check!
    void bksubPerm(Matrix A, Matrix B, Matrix P, Matrix result);
    //-------------------------LU factorization using Crout's Method--------------------------------//
    // factorizes the A matrix as the product of a unit upper triangular matrix U and a lower triangular matrix L
    uint8_t LU_Crout(Matrix A, Matrix L, Matrix U);
    //-------------------------LU factorization using Cormen's Method--------------------------------//
    // factorizes the A matrix as the product of a unit upper triangular matrix U and a lower triangular matrix L
    uint8_t LU_Cormen(Matrix A, Matrix L, Matrix U);
    //-----------------------LUP factorization using Cormen's Method------------------------------//
    // factorizes the A matrix as the product of a upper triangular matrix U and a unit lower triangular matrix L
    // returns the factor that has to be multiplied to the determinant of U in order to obtain the correct value
    int8_t LUP_Cormen(Matrix A, Matrix L, Matrix U, Matrix P);
    //-----------------------Linear system solver using LU factorization---------------------------//
    // solves the linear system A*X=B, where A is a n-by-n matrix and B an n-by-m matrix, giving the n-by-m matrix X
    void LinSolveLU(Matrix A, Matrix B, Matrix result);
    //----------------------Linear system solver using LUP factorization--------------------------//
    // solves the linear system A*X=B, where A is a n-by-n matrix and B an n-by-m matrix, giving the n-by-m matrix X
    void LinSolveLUP(Matrix A, Matrix B, Matrix result);
    //------------Linear system solver using Gauss elimination with partial pivoting---------------//
    // solves the linear system A*X=B, where A is a n-by-n matrix and B an n-by-m matrix, giving the n-by-m matrix X
    void LinSolveGauss(Matrix A, Matrix B, Matrix result);
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
    Matrix GaussNewton_Sens_Cal_9(Matrix Data, float k, Matrix X0, uint16_t nmax, float tol);
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
    Matrix GaussNewton_Sens_Cal_6(Matrix Data, float k, Matrix X0, uint16_t nmax, float tol);
    //------------------Quadratic form (sort of)----------------------//
    // returns matrix C=A*B*(~A)
    void QuadProd(Matrix A, Matrix B, Matrix result);
    
#ifdef __cplusplus
}
#endif

#endif /* NumMethodsF_h */
