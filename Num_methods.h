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
template <typename T, typename T2> Matrix<T> fwsub(Matrix<T> &A, Matrix<T2> &B){
    Matrix<T> result(A.columns(),B.columns());
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
}

//-------------------Backward substitution----------------------//
// assumes that the matrix A is already an upper triangular one. No check!
template <typename T, typename T2> Matrix<T> bksub(Matrix<T> &A, Matrix<T2> &B){
    int16_t ncolsA=A.columns();
    Matrix<T> result(ncolsA,B.columns());
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
}

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

template <typename T, typename T2> Matrix<T> GaussNewton_9(Matrix<T> &Data, Matrix<T2> &X0, uint16_t nmax, double tol){
    Matrix<T> result(X0);
    uint16_t nrows=Data.rows();
    uint16_t ncols=Data.columns();
    
    if ((nrows<9)||(ncols!=3))
        return Matrix<T>(0,0);
    
    Matrix<T> Jr(nrows,9);
    Matrix<T> res(nrows,1);
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
        Matrix<T> delta(Jr.pseudo_inv()*res);
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

template <typename T, typename T2> Matrix<T> GaussNewton_6(Matrix<T> &Data, Matrix<T2> &X0, uint16_t nmax, double tol){
    Matrix<T> result(X0);
    
    uint16_t nrows=Data.rows();
    uint16_t ncols=Data.columns();
    
    if ((nrows<6)||(ncols!=3))
        return Matrix<T>(0,0);
    
    Matrix<T> Jr(nrows,6);
    Matrix<T> res(nrows,1);
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
        Matrix<T> delta(Jr.pseudo_inv()*res);
        result-=delta;
        if (delta.norm()<tol)
            return result;
    }
    return result;
};

//------------Linear system solver using Gauss elimination with partial pivoting---------------//
// solves the linear system A*X=B, where A is a n-by-n matrix and B an n-by-m matrix, giving the n-by-m matrix X

template<typename T>  Matrix<T> linsolve(Matrix<T> &A, Matrix<T> &B) {
    int16_t pivrow=0;     // keeps track of current pivot row
    int16_t k;
    int16_t ncolsA=A.columns();
    int16_t ncolsB=B.columns();
    unsigned int i,j;      // k: overall index along diagonals; i: row index; j: col index
    double tmp;      // used for finding max value and making row swaps
    double tmp2;  // used to store abs when finding max value and to store coefficient value when eliminating values
    
    for (k = 0; k < (ncolsA-1); k++){
        
        // find pivot row, the row with biggest entry in current column
        tmp = fabs(A.get(k,k));
        pivrow=k;
        for (i = k+1; i < ncolsA; i++)
        {
            tmp2=fabs(A.get(i,k)); // 'Avoid using other functions inside abs()?'
            if (tmp2 >= tmp)
            {
                tmp = tmp2;
                pivrow = i;
            }
        }
        
        // check for singular Matrix
        if (A.get(pivrow,k) == 0.0)
        {
            return Matrix<T>(0,0);
        }
        
        // Execute pivot (row swap) if needed
        if (pivrow != k)
        {
            // swap row k of matrix A with pivrow
            for (j = 0; j < ncolsA; j++)
            {
                tmp = A.get(k, j);
                A.set(k, j) = A.get(pivrow,j);
                A.set(pivrow, j) = tmp;
            }
            // swap row k of matrix B with pivrow
            for (j = 0; j < ncolsB; j++)
            {
                tmp = B.get(k, j);
                B.set(k, j) = B.get(pivrow,j);
                B.set(pivrow, j) = tmp;
            }
        }
        
        // Row reduction
        tmp = 1.0/A.get(k,k);    // invert pivot element
        for (i=k+1;i<ncolsA;i++){ // along rows
            tmp2=A.get(i,k)*tmp;
            // Perform row reduction of A
            for (j=k+1;j<ncolsA;j++){ //along columns of A
                A.set(i,j)-=tmp2*A.get(k,j);
            }
            // Perform row reduction of B
            for (j=0;j<ncolsB;j++){ //along columns of B
                B.set(i,j)-=tmp2*B.get(k,j);
            }
        }
        
    }
    return(Matrix<T>(bksub(A,B)));
}


#endif
