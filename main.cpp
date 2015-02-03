//
//  main.cpp
//  test
//
//  Created by Andrea Vivani on 28/1/15.
//  Copyright (c) 2015 Andrea Vivani. All rights reserved.
//

#include <iostream>
#include "Matlib.h"

template <typename T> void mprint(MatrixX<T>& r) {
    std::cout << "\n" << (int) r.rows() << "x" << (int) r.columns() << "\n";
    for(int i=0; i<r.rows(); i++) {
        for (int j=0; j<r.columns(); j++)
            printf("%4.4f\t",(double) r(i,j));
        std::cout << "\n";
    }
}

/*class Test{
 public:
 Test (int a,int b){
     _a=a;
     _b=b;
 }
    int operator=(const std::initializer_list<int> &v){
        std::initializer_list<int>::iterator it;
        delete [] _ciccio;
        _ciccio=0;
        this->_ciccio= new int [v.size()];
        int ii=0;
        for (it=v.begin(); it!=v.end();it++)
        {
            _ciccio[ii]= *it;
            printf("%d\n",_ciccio[ii]);
            ii++;
        }
        return 1;
    }
    int* _ciccio;
    
private:
    int _a;
    int _b;
 };*/

int main(int argc, const char * argv[]) {
    // TEST GAUSS
     MatrixXd Data(10,3);
    Data={-0.009345794,	0.00952381,	1,
        -1,	-0.00952381,	0.044247788,
        1,	0.028571429,	0.008849558,
        -0.009345794,	1,	-0.115044248,
        0.028037383,	-1,	0.008849558,
        -0.775700935,	-0.638095238,	0.008849558,
        -0.570093458,	0.80952381,	-0.026548673,
        0.644859813,	0.771428571,	-0.061946903,
        0.775700935,	-0.619047619,	-0.115044248,
        0.981308411,	0.00952381,	-0.256637168};
    int8_t _data0[]={0,0,0,1,0,0,1,0,1};
     MatrixXs X0(9,1,_data0);
     int8_t _data01[]={0,0,0,1,1,1};
     MatrixXs X01(6,1,_data01);
     MatrixXd P=GaussNewton_9(Data,X0, 200,1e-9);
     MatrixXd P2=GaussNewton_6(Data, X01, 200, 1e-9);
     mprint(Data);
     mprint(P);
     mprint(P2);
    // TEST FWSUB
    double _dataMat1[]={  3,     0,     0,     0,
        5,     1,     0,    0,
        4,     5,     4,     0,
        5,     5,     2,     1};
    MatrixXd A1(4,4,_dataMat1);
    double _dataB1[]={ 0.276922984960890, 0.5, 0.046171390631154, 0.3, 0.097131781235848, 0.4, 0.823457828327293, 0.5};
    MatrixXd B1(4,2,_dataB1);
    MatrixXd S1=fwsub(A1,B1);
    mprint(S1);
    MatrixXd S3=A1*S1;
    mprint(S3);
    // TEST BKSUB
    double _dataMat2[]={3,     2,     1,     3,
        0,     2,     2,     4,
        0,     0,     2,     2,
        0,     0,     0,     3};
    MatrixXd A2(4,4,_dataMat2);
    MatrixXd S2=bksub(A2,B1);
    mprint(S2);
    MatrixXd S4=A2*S2;
    mprint(S4);
    
    // TEST LINSOLVE
    double _dataA[]={153,    0.3171,    0.3816,    0.4898,
        0.0462,    0.4358,    445,    0.4456,
        0.8235,    155,    0.7952,    0.6463,
        0.6948,    0.9745,    0.1869,    4456};
    double _dataB[]={0.7547,    0.1626,    0.3404,    0.2551,
        0.2760,    0.1190,    0.5853,    0.5060,
        0.6797,    0.4984,    0.2238,    0.6991,
        0.6551,    0.9597,    0.7513,   0.8909};
    MatrixXd AA(4,4,_dataA);
    MatrixXd BB(4,4,_dataB);
    MatrixXd Soli=!AA*BB;
    MatrixXd Solg=LinSolveGauss(AA, BB);
    MatrixXd Soll=LinSolveLU(AA, BB);
    MatrixXd Soll2=LinSolveLUP(AA, BB);
    mprint(Soli);
    mprint(Solg);
    mprint(Soll);
    mprint(Soll2);
    printf("\n%f\n",AA.det());

    
    
    return 0;
}