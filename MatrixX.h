//
//  MatrixX.h
//
//
//  Created by Andrea Vivani on 31/1/15.
//  Copyright (c) 2015 Andrea Vivani. All rights reserved.
//

#ifndef _MatrixX_H_
#define _MatrixX_H_
#include <stdint.h>

template <typename T>
class MatrixX {
public:
    //==================Assignment===================//
    MatrixX<T>(uint8_t m=0, uint8_t n=0, T* data=0);// constructor
    MatrixX<T>(const MatrixX<T> &rhs); // initializes one Matrix as equal to another, same datatype
    template <typename T2>  MatrixX<T>(const MatrixX<T2> &rhs); // initializes one Matrix as equal to another, different datatype
    inline T& set(uint8_t i, uint8_t j); // sets one single element
    inline T& operator()(uint8_t i, uint8_t j=0); //sets and returns one single element
    virtual ~MatrixX<T>(); // deconstructor
    MatrixX<T>& identity(); // sets the Matrix to an identity Matrix
    MatrixX<T>& zeros(); // fills the Matrix with zeros
    inline MatrixX<T>& operator=(const MatrixX<T> &rhs); // assignment
    MatrixX<T>& operator=(const std::initializer_list<T> &rhs); // assignment with list
    //==================Operations==================//
    template <typename T2> MatrixX<T> operator+(const MatrixX<T2> &rhs) const; // Matrix addition
    template <typename T2> MatrixX<T> operator+(T2 scalar) const; // Matrix and scalar addition
    template <typename T2> MatrixX<T>& operator+=(const MatrixX<T2> &rhs); // Matrix addition in place
    template <typename T2> MatrixX<T>& operator+=(T2 scalar); // Matrix scalar addition in place
    template <typename T2> MatrixX<T> operator-(const MatrixX<T2> &rhs) const; // Matrix subtraction
    template <typename T2> MatrixX<T> operator-(T2 scalar) const; // Matrix and scalar subtraction
    template <typename T2> MatrixX<T>& operator-=(const MatrixX<T2> &rhs); // Matrix subtraction in place
    template <typename T2> MatrixX<T>& operator-=(T2 scalar); // Matrix scalar subtraction in place
    template <typename T2> MatrixX<T> operator*(const MatrixX<T2> &rhs) const; // Matrix multiplication
    template <typename T2> MatrixX<T> operator*(T2 scalar) const; // Matrix scalar multiplication
    template <typename T2> MatrixX<T>& operator*=(const MatrixX<T2> &rhs);  // Matrix product in place
    template <typename T2> MatrixX<T>& operator*=(T2 scalar); // Matrix scalar multiplication in place
    template <typename T2> MatrixX<T> operator%(const MatrixX<T2> &rhs) const; //element-wise multiplication
    template <typename T2> MatrixX<T>& elw_multSelf(const MatrixX<T2> &rhs); //element-wise musltiplication in place
    template <typename T2> MatrixX<T> operator/(const MatrixX<T2> &rhs) const; //element-wise division
    template <typename T2> MatrixX<T>& elw_divSelf(const MatrixX<T2> &rhs); //element-wise musltiplication in place
    template <typename T2> MatrixX<T> cross(const MatrixX<T2> &rhs) const; //vector cross product
    MatrixX<T> operator-() const; // opposed Matrix
    MatrixX<T> operator!() const; // inversed Matrix
    MatrixX<T> inversed() const; //inversed Matrix using LU decomposition
    MatrixX<T> inversed_rob() const; //robust inversed Matrix using LUP decomposition
    MatrixX<T>& inverse(); //inverse in-place
    MatrixX<T> operator~() const; //transposed Matrix
    MatrixX<T> normalized(); //normalized Matrix
    MatrixX<T>& normalize(); //normalizes the Matrix in place
    MatrixX<T> pseudo_inv(); //Moore-Penrose pseudo inverse

    //==============Logical operations===============//
    template <typename T2> bool  operator==(const MatrixX<T2> &other) const; // true if two matrices are equal
    template <typename T2> bool  operator!=(const MatrixX<T2> &other) const; // true if two matrices are different
    
    //==================Matrix Data==================//
    inline T& get(uint8_t i, uint8_t j) const; //retrieves one single element
    T trace() const; // retrieves the trace
    inline uint8_t rows() const; // retrieves the number of rows
    inline uint8_t columns() const; // retrieves the number of columns
    inline uint8_t length() const; // retrieves the lenght of the Matrix
    double norm() const; //retrives the norm of the Matrix
    T sum() const; //retrieves the sum of all the elements
    T product()const; //retrieves the product of all the elements
    T det() const; //retrieves the determinant of the Matrix
    MatrixX<T> subMatrix(uint8_t row_top, uint8_t col_left, uint8_t row_bottom, uint8_t col_right) const; //returns a subMatrix of the original
    //===============Auxiliary Functions==============//
    template <typename T2> inline MatrixX<T>& copyData(const T2* data); //sets new data
    inline T& getData(uint8_t i) const; //retrieve a single data
    inline T* data() const; // retrieve all data
    
private:
    //-------Private Functions------//
    void release();
    void allocate();
    template <typename T2> MatrixX<T>& copyMatrix(const MatrixX<T2> &another);
    inline uint8_t index(uint8_t i, uint8_t j) const;
    //-------Private Variables-----//
    T* _data;
    uint8_t _nrows;
    uint8_t _ncols;
    bool _isAllocated=true;
};


typedef MatrixX<int8_t>                 MatrixXs;
typedef MatrixX<uint8_t>                MatrixXus;
typedef MatrixX<int16_t>                MatrixXi;
typedef MatrixX<uint16_t>               MatrixXui;
typedef MatrixX<int32_t>                MatrixXl;
typedef MatrixX<uint32_t>               MatrixXul;
typedef MatrixX<float>                  MatrixXf;
typedef MatrixX<double>                 MatrixXd;


#endif /* _MatrixX_H_ */
