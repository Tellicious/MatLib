//
//  MatrixX.tpp
//
//
//  Created by Andrea Vivani on 31/1/15.
//  Copyright (c) 2015 Andrea Vivani. All rights reserved.
//

#include "MatrixX.h"
#include <math.h>
#include <stdint.h>

//====================================Functions Prototypes=========================================//

template <typename T> bool LU_Crout(const MatrixX<T> &A, MatrixX<T> &L, MatrixX<T> &U);
template <typename T> bool LU_Cormen(const MatrixX<T> &A, MatrixX<T> &L, MatrixX<T> &U);
template <typename T> int8_t LUP_Cormen(const MatrixX<T> &A, MatrixX<T> &L, MatrixX<T> &U, MatrixXs &P);
template<typename T>  MatrixX<T> LinSolveLU(const MatrixX<T> &A, const MatrixX<T> &B);
template<typename T>  MatrixX<T> LinSolveLUP(const MatrixX<T> &A, const MatrixX<T> &B);
template<typename T>  MatrixX<T> LinSolveGauss(const MatrixX<T> &A, const MatrixX<T> &B);

//==========================================Assignment=============================================//
//-----------------------Constructor------------------------//
template<typename T>
MatrixX<T>::MatrixX(uint8_t m, uint8_t n, T* data) {
    this->_nrows = m;
    this->_ncols = n;
    this->_data = 0;
    allocate();
    if (data) {
        copyData(data);
    }
    else {
        zeros();
    }
}

//-----------Initialize as equal, same datatype-------------//
template<typename T>
MatrixX<T>::MatrixX(const MatrixX<T> &rhs) {
    _nrows = 0;
    _ncols = 0;
    _data = 0;
    copyMatrix(rhs);
}

//-----------Initialize as equal, diff datatype--------------//
template<typename T>
template <typename T2> MatrixX<T>::MatrixX(const MatrixX<T2> &rhs) {
    _nrows = 0;
    _ncols = 0;
    _data = 0;
    copyMatrix(rhs);
}

//--------------------Sets single element-------------------//
template<typename T>
T& MatrixX<T>::set(uint8_t i, uint8_t j) {
    return _data[(i * _ncols + j)];
}

//--------------Sets and returns single element-------------//
template<typename T>
T& MatrixX<T>::operator()(uint8_t i, uint8_t j){
    return set(i, j);
}

//----------------------Deconstructor-----------------------//
template<typename T>
MatrixX<T>::~MatrixX() {
    release();
}

//---------------------Identity Matrix----------------------//
template<typename T>
MatrixX<T>& MatrixX<T>::identity() {
    zeros();
    for(uint8_t i=0; i<_nrows; i++) {
        set(i,i) = 1.0;
    }
    return *this;
}

//----------------------Zeros Matrix-----------------------//
template<typename T>
MatrixX<T>& MatrixX<T>::zeros(){
    for (uint8_t i=0; i<_nrows * _ncols; i++) {
        this->_data[i] = 0.0;
    }
    return *this;
}

//-----------------------Assignment------------------------//
template<typename T>
MatrixX<T>& MatrixX<T>::operator=(const MatrixX<T> &rhs) {
    return copyMatrix(rhs);
}

//-----------------Assignment with list--------------------//
template<typename T>
MatrixX<T>& MatrixX<T>::operator=(const std::initializer_list<T> &rhs){
    //std::initializer_list<T>::iterator it;
    uint8_t ii=0;
    for (auto it : rhs){
        _data[ii]=it;
        ii++;
    }
    return *this;
}


//==========================================Operations=============================================//
//--------------------Matrix addition----------------------//
template<typename T>
template <typename T2> MatrixX<T> MatrixX<T>::operator+(const MatrixX<T2> &rhs) const{
    return MatrixX<T>(*this) += rhs;
}

//--------------------Scalar addition----------------------//
template<typename T>
template <typename T2> MatrixX<T> MatrixX<T>::operator+(T2 scalar) const{
    return MatrixX<T>(*this) += scalar;
}

//---------------Matrix addition in place-----------------//
template<typename T>
template <typename T2> MatrixX<T>& MatrixX<T>::operator+=(const MatrixX<T2> &rhs) {
    for (uint8_t i=0; i<_nrows*_ncols; i++)
        _data[i]+= rhs.getData(i);
    return *this;
}

//---------------Scalar addition in place-----------------//
template<typename T>
template <typename T2> MatrixX<T>& MatrixX<T>::operator+=(T2 scalar) {
    for (uint8_t i=0; i<_nrows*_ncols; i++)
        _data[i]+= scalar;
    return *this;
}
//------------------Matrix subtraction--------------------//
template<typename T>
template <typename T2>MatrixX<T> MatrixX<T>::operator-(const MatrixX<T2> &rhs) const{
    return MatrixX<T>(*this) -= rhs;
}
//------------------Scalar subtraction--------------------//
template<typename T>
template <typename T2> MatrixX<T> MatrixX<T>::operator-(T2 scalar) const{
    return MatrixX<T>(*this) -= scalar;
}

//-------------Matrix subtraction in place---------------//
template<typename T>
template <typename T2> MatrixX<T>& MatrixX<T>::operator-=(const MatrixX<T2> &rhs) {
    for (uint8_t i=0; i<_nrows*_ncols; i++)
        _data[i]-= rhs.getData(i);
    return *this;
}

//-------------Scalar subtraction in place---------------//
template<typename T>
template <typename T2> MatrixX<T>& MatrixX<T>::operator-=(T2 scalar) {
    for (uint8_t i=0; i<_nrows*_ncols; i++)
        _data[i]-= scalar;
    return *this;
}

//---------------Matrix multiplication------------------//
template<typename T>
template <typename T2> MatrixX<T> MatrixX<T>::operator*(const MatrixX<T2> &rhs) const {
    if (_ncols != rhs.rows()) {
        return MatrixX<T>(0, 0);
    }
    MatrixX<T> result(_nrows, rhs.columns());
    for (uint8_t i=0; i < _nrows; i++)
        for (uint8_t j=0; j < rhs.columns(); j++)
            for (uint8_t k=0; k < _ncols; k++)
                result(i, j) += get(i, k) * rhs.get(k, j);
    return result;
}

//---------------Scalar multiplication------------------//
template<typename T>
template <typename T2> MatrixX<T> MatrixX<T>::operator*(T2 scalar) const {
    return MatrixX<T>(*this)*=scalar;
}

//-----------Matrix multiplication in place-------------//
template<typename T>
template <typename T2> MatrixX<T>& MatrixX<T>::operator*=(const MatrixX<T2> &rhs){
    *this=*this*rhs;
    return *this;
}

//-----------Scalar multiplication in place-------------//
template<typename T>
template <typename T2> MatrixX<T>& MatrixX<T>::operator*=(T2 scalar){
    for (uint8_t i=0; i<_nrows * _ncols; i++)
        _data[i] *= scalar;
    return *this;
}

//-----------Element-wise multiplication ---------------//
template<typename T>
template <typename T2> MatrixX<T> MatrixX<T>::elw_mult(const MatrixX<T2> &rhs) const{
    // element-wise multiplication
    return MatrixX<T>(*this).elw_multSelf(rhs);
}

//--------Element-wise multiplication in place----------//
template<typename T>
template <typename T2> MatrixX<T>& MatrixX<T>::elw_multSelf(const MatrixX<T2> &rhs) {
    // element-wise multiplication with self-modification
    for (uint8_t i=0; i<_nrows*_ncols; i++)
        _data[i]*= rhs.getData(i);
    return *this;
}

//---------------Vector cross-product-------------------//
template<typename T>
template <typename T2> MatrixX<T> MatrixX<T>::cross(const MatrixX<T2> &rhs) const {
    if (((_nrows*_ncols)!=3) || ((rhs.columns()*rhs.rows())!=3)) //three element vectors only
        return MatrixX<T>(0,0); //empty
    MatrixX<T> result(3,1);
    result(0,0) = _data[1] * rhs.getData(2) - _data[2] * rhs.getData(1);
    result(1,0) = _data[2] * rhs.getData(0) - _data[0] * rhs.getData(2);
    result(2,0) = _data[0] * rhs.getData(1) - _data[1] * rhs.getData(0);
    return result;
}

//---------------------Opposed-------------------------//
template<typename T>
MatrixX<T> MatrixX<T>::operator-() const{
    return MatrixX<T>(*this) *= -1.0;
}

//---------------------Inverse-------------------------//
template<typename T>
MatrixX<T> MatrixX<T>::operator!() const{
    MatrixX<T> result(*this);
    result.inverse();
    return result;
}

//-----------------Inverse in place--------------------//
template<typename T>
MatrixX<T>& MatrixX<T>::inverse() {
    int16_t pivrow=0;     // keeps track of current pivot row
    int16_t k;
    unsigned int i,j;      // k: overall index along diagonal; i: row index; j: col index
    unsigned int* pivrows = new unsigned int[_ncols]; // keeps track of rows swaps to undo at end
    T tmp;      // used for finding max value and making column swaps
    T tmp_abs; // used to store abs when winding max value
    
    for (k = 0; k < _ncols; k++)
    {
        // find pivot row, the row with biggest entry in current column
        tmp = fabs(get(k,k));
        pivrow=k;
        for (i = k+1; i < _ncols; i++)
        {
            tmp_abs=fabs(get(i,k)); // 'Avoid using other functions inside abs()?'
            if (tmp_abs >= tmp)
            {
                tmp = tmp_abs;
                pivrow = i;
            }
        }
        
        // check for singular Matrix
        if (get(pivrow,k) == 0.0)
        {
            delete[] pivrows;
            this->release();
            return *this;
        }
        
        // Execute pivot (row swap) if needed
        if (pivrow != k)
        {
            // swap row k with pivrow
            for (j = 0; j < _ncols; j++)
            {
                tmp = get(k, j);
                set(k, j) = get(pivrow,j);
                set(pivrow, j) = tmp;
            }
        }
        pivrows[k] = pivrow;    // record row swap (even if no swap happened)
        
        tmp = 1.0/get(k,k);    // invert pivot element
        set(k,k) = 1.0;        // This element of input Matrix becomes result Matrix
        
        // Perform row reduction (divide every element by pivot)
        for (j = 0; j < _ncols; j++)
        {
            set(k,j) *= tmp;
        }
        
        // Now eliminate all other entries in this column
        for (i = 0; i < _ncols; i++)
        {
            if (i != k)
            {
                tmp = get(i,k);
                set(i,k) = 0.0;  // The other place where in Matrix becomes result mat
                for (j = 0; j < _ncols; j++)
                {
                    set(i,j) -= get(k,j) * tmp;
                }
            }
        }
    }
    
    // Done, now need to undo pivot row swaps by doing column swaps in reverse order
    for (k = _ncols-1; k >= 0; k--)
    {
        if (pivrows[k] != k)
        {
            for (i = 0; i < _ncols; i++)
            {
                tmp = get(i,k);
                set(i, k) = get(i,pivrows[k]);
                set(i,pivrows[k]) = tmp;
            }
        }
    }
    delete[] pivrows;
    return *this;
}

//-----------------Transposed--------------------//
template<typename T>
MatrixX<T> MatrixX<T>::operator~() const {
    MatrixX<T> result(_ncols,_nrows);
    for (uint8_t ii=0;ii<_nrows;ii++){
        for (uint8_t jj=0;jj<_ncols;jj++){
            result.set(jj,ii)=get(ii,jj);
        }
    }
    return result;
}

//-----------------Nomalized--------------------//
template<typename T>
MatrixX<T> MatrixX<T>::normalized() {
    return MatrixX<T>(*this).normalize();
}

//------------Nomalized in place---------------//
template<typename T>
MatrixX<T>& MatrixX<T>::normalize() {
    T k = norm();
    if (k > 0) {
        *this *= (1 / k);
    }
    return *this;
}

//-------Moore-Penrose pseudo inverse---------//
template<typename T>
MatrixX<T> MatrixX<T>::pseudo_inv() {
    MatrixX<T> tran=~(*this);
    return LinSolveLU((tran*(*this)), tran);
}

//=====================================Logical Operations=======================================//
//---------------True if equal----------------//
template<typename T>
template <typename T2> bool  MatrixX<T>::operator==(const MatrixX<T2> &other) const{
    if (_nrows != other.rows() || _ncols != other.columns())
        return false;
    
    for (uint8_t i=0; i<_nrows*_ncols; i++){
        if (_data[i] != other.getData(i)){
            return false;
        }
    }
    return true;
}

//-------------True if not equal--------------//
template<typename T>
template <typename T2> bool  MatrixX<T>::operator!=(const MatrixX<T2> &other) const{
    return !((*this) == other);
}

//=======================================Matrix Data=========================================//
//---------Returns one single element---------//
template<typename T>
const T& MatrixX<T>::get(uint8_t i, uint8_t j) const{
    return _data[(i * _ncols + j)];
}

//-------------Returns the trace-------------//
template<typename T>
T MatrixX<T>::trace() const {
    T result = 0.0;
    for (uint8_t i=0; i<_nrows && i<_ncols; i++) {
        result += get(i, i);
    }
    return result;
}

//---------Returns the number of rows--------//
template<typename T>
uint8_t MatrixX<T>::rows() const{
    return _nrows;
}

//-------Returns the number of columns-------//
template<typename T>
uint8_t MatrixX<T>::columns() const{
    return _ncols;
}

//------------Returns the length------------//
template<typename T>
uint8_t MatrixX<T>::length() const{
    return fmax(_nrows, _ncols);
}

//-------------Returns the norm--------------//
template<typename T>
double MatrixX<T>::norm() const {
    double result = 0.0;
    for (uint8_t i=0; i<_nrows*_ncols; i++){
        result += _data[i]*_data[i];
    }
    result = sqrt(result);
    return result;
}

//--------------Returns the sum--------------//
template<typename T>
T MatrixX<T>::sum() const {
    T result = 0.0;
    for (uint8_t i=0; i<_nrows*_ncols; i++)
        result += _data[i];
    return result;
}

//------------Returns the product------------//
template<typename T>
T MatrixX<T>::product() const {
    T result = 1.0;
    for (uint8_t i=0; i<_nrows*_ncols; i++)
        result *= _data[i];
    return result;
}
//-----------Returns the determinant----------//
template<typename T>
T MatrixX<T>::det() const{
    MatrixX<T> U(_nrows,_nrows);
    MatrixX<T> L(_nrows,_nrows);
    int16_t ii;
    T determinant=1;
    if(LU_Cormen(*this, L, U)){
        for (ii=0;ii<_nrows;ii++){
            determinant*=U.get(ii,ii);
        }
        return determinant;
    }
    else{
        MatrixXs P(_nrows,1);
        int8_t det_f=LUP_Cormen(*this,L,U,P);
        if(det_f){
            for (ii=0;ii<_nrows;ii++){
                determinant*=U.get(ii,ii);
            }
            return determinant*det_f;
        }
        else {
            return 0.0;
        }
    }
}

//------------Returns a subMatrix------------//
template<typename T>
MatrixX<T> MatrixX<T>::subMatrix(uint8_t row_top, uint8_t col_left, uint8_t row_bottom, uint8_t col_right) const {
    uint8_t rows = row_bottom-row_top+1;
    uint8_t cols = col_right-col_left+1;
    MatrixX<T> result(rows, cols);
    
    for(uint8_t i=0;i<rows; i++)
        for(uint8_t j=0; j<cols; j++)
            result(i,j) = get(row_top+i, col_left+j);
    return result;
    
}

//==========================================Auxiliary=============================================//
//------------Copies data of one Matrix to another-------------//
template<typename T>
template <typename T2> MatrixX<T>& MatrixX<T>::copyData(const T2* data) {
    for (uint8_t i=0; i<_nrows*_ncols; i++)
        this->_data[i] = (T) data[i];
    return *this;
}

//------------Returns one single value of data-----------------//
template<typename T>
T& MatrixX<T>::getData(uint8_t i) const{
    return _data[i];
}

//-------------Returns the whole data vector-------------------//
template<typename T>
T* MatrixX<T>::data() const{
    return this->_data;
}

//===========================================Private==============================================//
//----------------------Releases memory------------------------//
template<typename T>
void MatrixX<T>::release() {
    if (_data && _isAllocated) {
        delete[] _data;
        _data = 0;
        _isAllocated = false;
    }
}

//----------------------Allocates memory-----------------------//
template<typename T>
void MatrixX<T>::allocate() {
    release();
    if (_nrows && _ncols)
        this->_data = new T[_nrows * _ncols];
    else
        this->_data = 0;
    _isAllocated = true;
}

//----------------Copies one Matrix to another-----------------//
template<typename T>
template <typename T2> MatrixX<T>& MatrixX<T>::copyMatrix(const MatrixX<T2> &another) {
    if ((*this) != another) {
        if (_nrows * _ncols != another.rows() * another.columns()) {
            _nrows = another.rows();
            _ncols = another.columns();
            allocate();
        }
        else{
            _nrows = another.rows();
            _ncols = another.columns();
        }
        copyData(another.data());
    }
    return *this;
}

//---------------Data index from row and column----------------//
template<typename T>
uint8_t MatrixX<T>::index(uint8_t i, uint8_t j) const{
    return i * _ncols + j;
}

