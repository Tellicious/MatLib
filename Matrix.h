//
//  Matrix.h
//
//
//  Created by Andrea Vivani on 31/1/15.
//  Copyright (c) 2015 Andrea Vivani. All rights reserved.
//

#ifndef Matrix_H_
#define Matrix_H_
#include <stdint.h>

template <typename T>
class Matrix {
public:
    //==================Assignment===================//
    Matrix<T>(uint8_t m=0, uint8_t n=0, T* data=0);// constructor
    Matrix<T>(const Matrix<T> &rhs); // initializes one Matrix as equal to another, same datatype
    template <typename T2>  Matrix<T>(const Matrix<T2> &rhs); // initializes one Matrix as equal to another, different datatype
    T& set(uint8_t i, uint8_t j); // sets one single element
    T& operator()(uint8_t i, uint8_t j=0); //sets and returns one single element
    virtual ~Matrix<T>(); // deconstructor
    Matrix<T>& identity(); // sets the Matrix to an identity Matrix
    Matrix<T>& zeros(); // fills the Matrix with zeros
    Matrix<T>& operator=(const Matrix<T> &rhs); // assignment
    Matrix<T>& operator=(const std::initializer_list<T> &rhs); // assignment with list
    //==================Operations==================//
    template <typename T2> Matrix<T> operator+(const Matrix<T2> &rhs) const; // Matrix addition
    template <typename T2> Matrix<T> operator+(T2 scalar) const; // Matrix and scalar addition
    template <typename T2> Matrix<T>& operator+=(const Matrix<T2> &rhs); // Matrix addition in place
    template <typename T2> Matrix<T>& operator+=(T2 scalar); // Matrix scalar addition in place
    template <typename T2> Matrix<T> operator-(const Matrix<T2> &rhs) const; // Matrix subtraction
    template <typename T2> Matrix<T> operator-(T2 scalar) const; // Matrix and scalar subtraction
    template <typename T2> Matrix<T>& operator-=(const Matrix<T2> &rhs); // Matrix subtraction in place
    template <typename T2> Matrix<T>& operator-=(T2 scalar); // Matrix scalar subtraction in place
    template <typename T2> Matrix<T> operator*(const Matrix<T2> &rhs) const; // Matrix multiplication
    template <typename T2> Matrix<T> operator*(T2 scalar) const; // Matrix scalar multiplication
    template <typename T2> Matrix<T>& operator*=(const Matrix<T2> &rhs);  // Matrix product in place
    template <typename T2> Matrix<T>& operator*=(T2 scalar); // Matrix scalar multiplication in place
    template <typename T2> Matrix<T> elw_mult(const Matrix<T2> &rhs) const; //element-wise multiplication
    template <typename T2> Matrix<T>& elw_multSelf(const Matrix<T2> &rhs); //element-wise musltiplication in place
    template <typename T2> Matrix<T> cross(const Matrix<T2> &rhs) const; //vector cross product
    Matrix<T> operator-() const; // opposed Matrix
    Matrix<T> operator!() const; // inverse Matrix
    Matrix<T>& inverse(); //inverse in-place
    Matrix<T> operator~() const; //transposed Matrix
    Matrix<T> normalized(); //normalized Matrix
    Matrix<T>& normalize(); //normalizes the Matrix in place
    Matrix<T> pseudo_inv(); //Moore-Penrose pseudo inverse

    //==============Logical operations===============//
    template <typename T2> bool  operator==(const Matrix<T2> &other) const; // true if two matrices are equal
    template <typename T2> bool  operator!=(const Matrix<T2> &other) const; // true if two matrices are different
    
    //==================Matrix Data==================//
    const T& get(uint8_t i, uint8_t j) const; //retrieves one single element
    T trace() const; // retrieves the trace
    uint8_t rows() const; // retrieves the number of rows
    uint8_t columns() const; // retrieves the number of columns
    uint8_t length() const; // retrieves the lenght of the Matrix
    double norm() const; //retrives the norm of the Matrix
    T sum() const; //retrieves the sum of all the elements
    T product()const; //retrieves the product of all the elements
    //T det() const; //retrieves the determinant of the Matrix
    Matrix<T> subMatrix(uint8_t row_top, uint8_t col_left, uint8_t row_bottom, uint8_t col_right) const; //returns a subMatrix of the original
    //===============Auxiliary Functions==============//
    template <typename T2> Matrix<T>& copyData(const T2* data); //sets new data
    T& getData(uint8_t i) const; //retrieve a single data
    T* data() const; // retrieve all data
    
private:
    //-------Private Functions------//
    void release();
    void allocate();
    template <typename T2> Matrix<T>& copyMatrix(const Matrix<T2> &another);
    uint8_t index(uint8_t i, uint8_t j) const;
    //-------Private Variables-----//
    T* _data;
    uint8_t _nrows;
    uint8_t _ncols;
    bool _isAllocated=true;
};


typedef Matrix<int8_t>                 MatrixXs;
typedef Matrix<uint8_t>                MatrixXus;
typedef Matrix<int16_t>                MatrixXi;
typedef Matrix<uint16_t>               MatrixXui;
typedef Matrix<int32_t>                MatrixXl;
typedef Matrix<uint32_t>               MatrixXul;
typedef Matrix<float>                  MatrixXf;
typedef Matrix<double>                 MatrixXd;


#endif /* Matrix2_H_ */
