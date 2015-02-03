//
//  Num_methods.h
//
//
//  Created by Andrea Vivani on 31/1/15.
//  Copyright (c) 2015 Andrea Vivani. All rights reserved.
//

#ifndef test_Num_methods_h
#define test_Num_methods_h
#include <math.h>
#include <stdint.h>

//-------------------Forward substitution----------------------//
// assumes that the matrix A is already a lower triangular one. No check!

template <typename T, typename T2> MatrixX<T> fwsub(const MatrixX<T> &A, const MatrixX<T2> &B){
    MatrixX<T> result(A.columns(),B.columns());
    for (int k=0;k<B.columns();k++){
        result.set(0,k) = B.get(0,k) / A.get(0,0);
        for (int i=1;i<A.rows();i++){
            T tmp=0.0;
            for (int j=0;j<i;j++){
                tmp+=A.get(i,j) * (result.get(j,k));
            }
            result.set(i, k)= (B.get(i,k) - tmp) / A.get(i,i);
        }
    }
    return result;
};

//---------------Forward substitution with permutation-------------------//
// assumes that the matrix A is already a lower triangular one. No check!

template <typename T, typename T2> MatrixX<T> fwsub_P( const MatrixX<T> &A, const MatrixX<T2> &B, MatrixXs &P){
    MatrixX<T> result(A.columns(),B.columns());
    for (int k=0;k<B.columns();k++){
        result.set(0,k) = B.get(P.get(0,0),k) / A.get(0,0);
        for (int i=1;i<A.rows();i++){
            T tmp=0.0;
            for (int j=0;j<i;j++){
                tmp+=A.get(i,j) * (result.get(j,k));
            }
            result.set(i, k)= (B.get(P.get(i,0),k) - tmp) / A.get(i,i);
        }
    }
    return result;
};

//-------------------Backward substitution----------------------//
// assumes that the matrix A is already an upper triangular one. No check!

template <typename T, typename T2> MatrixX<T> bksub(const MatrixX<T> &A, const MatrixX<T2> &B){
    int16_t ncolsA=A.columns();
    MatrixX<T> result(ncolsA,B.columns());
    for (int k=0;k<B.columns();k++){
        result.set(ncolsA-1,k) = B.get(ncolsA-1,k) / A.get(ncolsA-1,ncolsA-1);
        for (int i=A.rows()-2;i>=0;i--){
            T tmp=0.0;
            for (int j=ncolsA-1;j>i;j--){
                tmp+=A.get(i,j) * (result.get(j,k));
            }
            result.set(i, k)= (B.get(i,k) - tmp) / A.get(i,i);
        }
    }
    return result;
};

//--------------Backward substitution with permutation-----------------//
// assumes that the matrix A is already an upper triangular one. No check!

template <typename T, typename T2> MatrixX<T> bksub_P(const MatrixX<T> &A, const MatrixX<T2> &B, MatrixXs &P){
    int16_t ncolsA=A.columns();
    MatrixX<T> result(ncolsA,B.columns());
    for (int k=0;k<B.columns();k++){
        result.set(ncolsA-1,k) = B.get(P.get(ncolsA-1,0),k) / A.get(ncolsA-1,ncolsA-1);
        for (int i=A.rows()-2;i>=0;i--){
            T tmp=0.0;
            for (int j=ncolsA-1;j>i;j--){
                tmp+=A.get(i,j) * (result.get(j,k));
            }
            result.set(i, k)= (B.get(P.get(i,0),k) - tmp) / A.get(i,i);
        }
    }
    return result;
};

//-------------------------LU factorization using Crout's Method--------------------------------//
// factorizes the A matrix as the product of a unit upper triangular matrix U and a lower triangular matrix L

template <typename T> bool LU_Crout(const MatrixX<T> &A, MatrixX<T> &L, MatrixX<T> &U){
    MatrixX<T> A_tmp(A);
    int16_t ii, jj, kk;
    T sum = 0;
    int16_t nrowsA=A.rows();
    U.identity();
    L.zeros();
    for (jj = 0; jj < nrowsA; jj++) {
        for (ii = jj; ii < nrowsA; ii++) {
            sum = 0;
            for (kk = 0; kk < jj; kk++) {
                sum += L.get(ii,kk) * U(kk,jj);
            }
            L.set(ii,jj) = A.get(ii,jj) - sum;
        }
        
        for (ii = jj; ii < nrowsA; ii++) {
            sum = 0;
            for(kk = 0; kk < jj; kk++) {
                sum += L.get(jj,kk) * U.get(kk,ii);
            }
            if (L.get(jj,jj) == 0) {
                return false;
            }
            U.set(jj,ii) = (A.get(jj,ii) - sum) / L.get(jj,jj);
        }
    }
    return true;
};

//-------------------------LU factorization using Cormen's Method--------------------------------//
// factorizes the A matrix as the product of a unit upper triangular matrix U and a lower triangular matrix L
template <typename T> bool LU_Cormen(const MatrixX<T> &A, MatrixX<T> &L, MatrixX<T> &U){
    MatrixX<T> A_tmp(A);
    int16_t i, j, k;
    T tmp;
    int16_t nrowsA=A.rows();
    L.identity();
    U.zeros();
    
    for (k=0; k<nrowsA; k++) {
        U.set(k,k)=A_tmp.get(k,k);
        if (A_tmp.get(k,k)==0){
            return false;
        }
        tmp=1.0/U.get(k,k);
        for (i=k+1; i<nrowsA; i++) {
            L.set(i,k)=A_tmp.get(i,k)*tmp;
            U.set(k,i)=A_tmp.get(k,i);
        }
        for (i=k+1; i<nrowsA; i++) {
            for (j=k+1; j<nrowsA; j++) {
                A_tmp.set(i,j)-=L.get(i,k)*U.get(k,j);
            }
        }
    }
    return true;
};

//-----------------------LUP factorization using Cormen's Method------------------------------//
// factorizes the A matrix as the product of a upper triangular matrix U and a unit lower triangular matrix L
// returns the factor that has to be multiplied to the determinant of U in order to obtain the correct value

template <typename T> int8_t LUP_Cormen(const MatrixX<T> &A, MatrixX<T> &L, MatrixX<T> &U, MatrixXs &P){
    MatrixX<T> A_tmp(A);
    int16_t nrowsA=A_tmp.rows();
    int16_t i,j,k;
    T tmp,tmp2;
    int8_t d_mult=1; // determinant multiplying factor
    // initialization
    for (i = 0; i < nrowsA; i++){
        P(i,0) = i;
    }
    U.zeros();
    L.identity();
    
    // outer loop over diagonal pivots
    for (k = 0; k < nrowsA-1; k++) {
        
        // inner loop to find the largest pivot
        int16_t pivrow = k;
        tmp=fabs(A_tmp.get(k,k));
        for (i = k+1 ; i < nrowsA; i++){
            tmp2=fabs(A_tmp.get(i,k));
            if (tmp2 > tmp){
                tmp=tmp2;
                pivrow = i;
            }
        }
        
        // check for singularity
        if (A_tmp.get(pivrow,k)==0) {
            return 0;
        }
        
        // swap rows
        if (pivrow != k) {
            P.set(k,0)=pivrow;
            P.set(pivrow,0)=k;
            d_mult*=-1;
            
            for (j = 0; j < nrowsA; j++){
                tmp=A_tmp.get(k,j);
                A_tmp.set(k,j)=A_tmp.get(pivrow,j);
                A_tmp.set(pivrow,j)=tmp;
            }
        }
        tmp=1.0/A_tmp.get(k,k);
        // Gaussian elimination
        for (i = k + 1; i < nrowsA; i++) { // iterate down rows
            A_tmp.set(i,k) *= tmp;
            for (j = k + 1; j < nrowsA; j++){ // iterate across rows
                A_tmp.set(i,j) -= A_tmp.get(i,k) * A_tmp.get(k,j);
            }
        }
    }
    for (k=0;k<nrowsA;k++){
        U.set(k,k)=A_tmp.get(k,k);
        for (j=k+1;j<nrowsA;j++){
            L.set(j,k)=A_tmp.get(j,k);
            U.set(k,j)=A_tmp.get(k,j);
        }
    }
    return d_mult;
};

//-----------------------Linear system solver using LU factorization---------------------------//
// solves the linear system A*X=B, where A is a n-by-n matrix and B an n-by-m matrix, giving the n-by-m matrix X

template<typename T>  MatrixX<T> LinSolveLU(const MatrixX<T> &A, const MatrixX<T> &B) {
    MatrixX<T> L(A.rows(),A.columns());
    MatrixX<T> U(L);
    LU_Crout(A, L, U);
    MatrixX<T> y=fwsub(L,B);
    return MatrixX<T>(bksub(U,y));
};

//----------------------Linear system solver using LUP factorization--------------------------//
// solves the linear system A*X=B, where A is a n-by-n matrix and B an n-by-m matrix, giving the n-by-m matrix X

template<typename T>  MatrixX<T> LinSolveLUP(const MatrixX<T> &A, const MatrixX<T> &B) {
    MatrixX<T> L(A.rows(),A.columns());
    MatrixX<T> U(L);
    MatrixXs P(A.rows(),1);
    LUP_Cormen(A, L, U, P);
    MatrixX<T> y=fwsub_P(L,B,P);
    return MatrixX<T>(bksub(U,y));
};

//------------Linear system solver using Gauss elimination with partial pivoting---------------//
// solves the linear system A*X=B, where A is a n-by-n matrix and B an n-by-m matrix, giving the n-by-m matrix X

template<typename T>  MatrixX<T> LinSolveGauss(const MatrixX<T> &A, const MatrixX<T> &B) {
    MatrixX<T> A_tmp=A;
    MatrixX<T> B_tmp=B;
    int16_t pivrow=0;     // keeps track of current pivot row
    int16_t ncolsA=A_tmp.columns();
    int16_t ncolsB=B_tmp.columns();
    int16_t k,i,j;      // k: overall index along diagonals; i: row index; j: col index
    double tmp;      // used for finding max value and making row swaps
    double tmp2;  // used to store abs when finding max value and to store coefficient value when eliminating values
    
    for (k = 0; k < (ncolsA-1); k++){
        
        // find pivot row, the row with biggest entry in current column
        tmp = fabs(A_tmp.get(k,k));
        pivrow=k;
        for (i = k+1; i < ncolsA; i++)
        {
            tmp2=fabs(A_tmp.get(i,k)); // 'Avoid using other functions inside abs()?'
            if (tmp2 >= tmp)
            {
                tmp = tmp2;
                pivrow = i;
            }
        }
        
        // check for singular Matrix
        if (A_tmp.get(pivrow,k) == 0.0)
        {
            return MatrixX<T>(0,0);
        }
        
        // Execute pivot (row swap) if needed
        if (pivrow != k)
        {
            // swap row k of matrix A with pivrow
            for (j = k; j < ncolsA; j++)
            {
                tmp = A_tmp.get(k, j);
                A_tmp.set(k, j) = A_tmp.get(pivrow,j);
                A_tmp.set(pivrow, j) = tmp;
            }
            // swap row k of matrix B with pivrow
            for (j = 0; j < ncolsB; j++)
            {
                tmp = B_tmp.get(k, j);
                B_tmp.set(k, j) = B_tmp.get(pivrow,j);
                B_tmp.set(pivrow, j) = tmp;
            }
        }
        
        // Row reduction
        tmp = 1.0/A_tmp.get(k,k);    // invert pivot element
        for (i=k+1;i<ncolsA;i++){ // along rows
            tmp2=A_tmp.get(i,k)*tmp;
            // Perform row reduction of A
            for (j=k+1;j<ncolsA;j++){ //along columns of A
                A_tmp.set(i,j)-=tmp2*A_tmp.get(k,j);
            }
            // Perform row reduction of B
            for (j=0;j<ncolsB;j++){ //along columns of B
                B_tmp.set(i,j)-=tmp2*B_tmp.get(k,j);
            }
        }
        
    }
    return(MatrixX<T>(bksub(A_tmp,B_tmp)));
};

//------------Gauss-Newton Method with 9 parameters---------------//
// approximates Data to a sphere by calculating 6 gains (s) and 3 biases (b), useful to calibrate some sensors (meas_sphere=S*(meas-B) with S symmetric)
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

template <typename T, typename T2> MatrixX<T> GaussNewton_9(const MatrixX<T> &Data, MatrixX<T2> &X0, uint16_t nmax, double tol){
    MatrixX<T> result(X0);
    uint16_t nrows=Data.rows();
    uint16_t ncols=Data.columns();
    
    if ((nrows<9)||(ncols!=3))
        return MatrixX<T>(0,0);
    
    MatrixX<T> Jr(nrows,9);
    MatrixX<T> res(nrows,1);
    for (int n_iter=1;n_iter<nmax;n_iter++){
        for (int jj=0;jj<nrows;jj++){
            T d1=Data.get(jj,0) - result.get(0,0);
            T d2=Data.get(jj,1) - result.get(1,0);
            T d3=Data.get(jj,2) - result.get(2,0);
            T rx1= -2*(result.get(3,0)*d1 + result.get(4,0)*d2 + result.get(5,0)*d3);
            T rx2= -2*(result.get(4,0)*d1 + result.get(6,0)*d2 + result.get(7,0)*d3);
            T rx3= -2*(result.get(5,0)*d1 + result.get(7,0)*d2 + result.get(8,0)*d3);
            Jr.set(jj,0)=result.get(3,0)*rx1+result.get(4,0)*rx2+result.get(5,0)*rx3;
            Jr.set(jj,1)=result.get(4,0)*rx1+result.get(6,0)*rx2+result.get(7,0)*rx3;
            Jr.set(jj,2)=result.get(5,0)*rx1+result.get(7,0)*rx2+result.get(8,0)*rx3;
            Jr.set(jj,3)=-d1*rx1;
            Jr.set(jj,4)=-d2*rx1-d1*rx2;
            Jr.set(jj,5)=-d3*rx1-d1*rx3;
            Jr.set(jj,6)=-d2*rx2;
            Jr.set(jj,7)=-d3*rx2-d2*rx3;
            Jr.set(jj,8)=-d3*rx3;
            T t1=result.get(3,0)*d1 + result.get(4,0)*d2 + result.get(5,0)*d3;
            T t2=result.get(4,0)*d1 + result.get(6,0)*d2 + result.get(7,0)*d3;
            T t3=result.get(5,0)*d1 + result.get(7,0)*d2 + result.get(8,0)*d3;
            res.set(jj,0)= t1*t1+t2*t2+t3*t3-1;
        }
        // pseudo-inverse*res (inv(Jr'*Jr)*Jr'*res)
        MatrixX<T> delta(Jr.pseudo_inv()*res);
        result-=delta;
        if (delta.norm()<tol)
            return result;
    }
    return result;
};


//------------Gauss-Newton Method with 6 parameters---------------//
// approximates Data to a sphere by calculating 3 gains (s) and 3 biases (b), useful to calibrate some sensors (meas_sphere=S*(meas-B) with S diagonal)
// Data has n>=6 rows corresponding to the number of measures and 3 columns corresponding to the 3 axis
// X0 is the starting guess vector (usually [0 0 0 1 1 1]), nmax the maximum number of iterations (200 is generally fine, even if it usually converges within 10 iterations), and tol the stopping tolerance (1e-6 is usually more than fine)
/*b1=out(0,0);
 b2=out(1,0);
 b3=out(2,0);
 s11=out(3,0);
 s22=out(4,0);
 s33=out(5,0);*/

template <typename T, typename T2> MatrixX<T> GaussNewton_6(const MatrixX<T> &Data, MatrixX<T2> &X0, uint16_t nmax, double tol){
    MatrixX<T> result(X0);
    
    uint16_t nrows=Data.rows();
    uint16_t ncols=Data.columns();
    
    if ((nrows<6)||(ncols!=3))
        return MatrixX<T>(0,0);
    
    MatrixX<T> Jr(nrows,6);
    MatrixX<T> res(nrows,1);
    for (int n_iter=1;n_iter<nmax;n_iter++){
        for (int jj=0;jj<nrows;jj++){
            T d1=Data.get(jj,0) - result.get(0,0);
            T d2=Data.get(jj,1) - result.get(1,0);
            T d3=Data.get(jj,2) - result.get(2,0);
            Jr.set(jj,0)=-2*d1*result.get(3,0)*result.get(3,0);
            Jr.set(jj,1)=-2*d2*result.get(4,0)*result.get(4,0);
            Jr.set(jj,2)=-2*d3*result.get(5,0)*result.get(5,0);
            Jr.set(jj,3)=2*result.get(3,0)*d1*d1;
            Jr.set(jj,4)=2*result.get(4,0)*d2*d2;
            Jr.set(jj,5)=2*result.get(5,0)*d3*d3;
            T t1=result.get(3,0)*d1;
            T t2=result.get(4,0)*d2;
            T t3=result.get(5,0)*d3;
            res.set(jj,0)= t1*t1+t2*t2+t3*t3-1;
        }
        // pseudo-inverse*res (inv(Jr'*Jr)*Jr'*res)
        MatrixX<T> delta(Jr.pseudo_inv()*res);
        result-=delta;
        if (delta.norm()<tol)
            return result;
    }
    return result;
};

#endif
