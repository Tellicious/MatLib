//
//  main.cpp
//  test
//
//  Created by Andrea Vivani on 28/1/15.
//  Copyright (c) 2015 Andrea Vivani. All rights reserved.
//

#include <iostream>
#include <chrono>
#include "Matlib.h"

template <typename T> void mprint(MatrixX<T>& r) {
    std::cout << "\n" << (int) r.rows() << "x" << (int) r.columns() << "\n";
    for(int i=0; i<r.rows(); i++) {
        for (int j=0; j<r.columns(); j++)
            printf("%4.4f\t",(double) r(i,j));
        std::cout << "\n";
    }
}

/*MatrixXd qqq(MatrixXd A, MatrixXd B){
    MatrixXd result(A.rows(),A.rows());
    //for (int ii=0; ii<A.rows(); ii++) {
        for (int jj=0; jj<A.rows(); jj++) {
            for (int i=0; i<A.columns(); i++) {
                double ttt=0;
                for (int j=0;j<A.columns();j++){
                    ttt+=A.set(jj,j)*B.get(i,j);
                }
                for (int ii=0; ii<A.rows(); ii++) {
                result.set(ii,jj)+=A.get(ii,i)*ttt;
            }
        }
    }
    return result;
}*/

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
    std::cout<<"Test Gauss\n";
    MatrixXd Data(20,3);
    /*Data={-0.009345794,	0.00952381,	1,
     -1,	-0.00952381,	0.044247788,
     1,	0.028571429,	0.008849558,
     -0.009345794,	1,	-0.115044248,
     0.028037383,	-1,	0.008849558,
     -0.775700935,	-0.638095238,	0.008849558,
     -0.570093458,	0.80952381,	-0.026548673,
     0.644859813,	0.771428571,	-0.061946903,
     0.775700935,	-0.619047619,	-0.115044248,
     0.981308411,	0.00952381,	-0.256637168};*/
    
    Data={0.6568, -0.6316, -9.7223,
        -0.1856, -0.1258, -9.4601,
        -0.1487, -1.9148, -9.1830,
        -1.2664, 1.3523, -2.5576,
        -3.0869, -0.9538, 9.3472,
        -1.9048, -2.8481, 9.4088,
        -0.7129, -0.2773, 9.5633,
        7.9806, 2.3157, 2.7211,
        1.7272, -4.6918, 0.5029,
        -7.1359, -6.8862, 0.5362,
        -7.1575, -6.8544, 0.6613,
        -7.5065, -1.7227, 0.7005,
        -6.7207, 7.4355, 0.5626,
        -4.1161, 4.1446, -4.4461,
        -0.1862, -0.1781, -9.8555,
        -0.1886, -0.1756, -9.8548,
        -0.1913, -0.1744, -9.8540,
        -0.1907, -0.1738, -9.8542,
        -0.1922, -0.1707, -9.8523,
        -0.1903, -0.1735, -9.8529};
    double _data0[]={0,0,0,1,0,0,1,0,1};
    MatrixXd X0(9,1,_data0);
    double _data01[]={0,0,0,1,1,1};
    MatrixXd X01(6,1,_data01);
    MatrixXd P=GaussNewton_9(Data,X0, 200,1e-9);
    MatrixXd P2=GaussNewton_6(Data, X01, 200, 1e-9);
    std::cout<<"\nData\n";
    mprint(Data);
    std::cout<<"\n9 Parametri:\n";
    mprint(P);
    std::cout<<"\n6 Parametri:\n";
    mprint(P2);
    // TEST LINSOLVE
    std::cout<<"\nTest LinSolve\n";
    double _dataA[]={0.5432,    0.3171,    0.3816,    0.4898,
        0.0462,    0.4358,    0.6651,    0.4456,
        0.8235,    0.1324,    0.7952,    0.6463,
        0.6948,    0.9745,    0.1869,    0.4456};
    double _dataB[]={0.7547,    0.1626,    0.3404,    0.2551,
        0.2760,    0.1190,    0.5853,    0.5060,
        0.6797,    0.4984,    0.2238,    0.6991,
        0.6551,    0.9597,    0.7513,   0.8909};
    MatrixXd AA(4,4,_dataA);
    MatrixXd BB(4,4,_dataB);
    auto t0 = std::chrono::high_resolution_clock::now();
    MatrixXd Soli=!AA*BB;
    auto t1 = std::chrono::high_resolution_clock::now();
    MatrixXd Solg=LinSolveGauss(AA, BB);
    auto t2 = std::chrono::high_resolution_clock::now();
    MatrixXd Soll=LinSolveLU(AA, BB);
    auto t3 = std::chrono::high_resolution_clock::now();
    MatrixXd Soll2=LinSolveLUP(AA, BB);
    auto t4 = std::chrono::high_resolution_clock::now();
    std::cout <<"Inversa: "<< std::chrono::duration_cast<std::chrono::microseconds>(t1-t0).count()<<"\n";
    std::cout <<"Gauss: "<< std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count()<<"\n";
    std::cout <<"LU: "<< std::chrono::duration_cast<std::chrono::microseconds>(t3-t2).count()<<"\n";
    std::cout <<"LUP: "<< std::chrono::duration_cast<std::chrono::microseconds>(t4-t3).count()<<"\n";
    std::cout <<"\nSoluz. Inversa:\n";
    mprint(Soli);
    std::cout <<"\nSoluz. Gauss:\n";
    mprint(Solg);
    std::cout <<"\nSoluz. LU:\n";
    mprint(Soll);
    std::cout <<"\nSoluz. LUP:\n";
    mprint(Soll2);
    // TEST INVERSA
    std::cout << "\nTest Inversa"<<"\n";
    auto t5 = std::chrono::high_resolution_clock::now();
    MatrixXd invA1=!AA;
    auto t6 = std::chrono::high_resolution_clock::now();
    MatrixXd invA2=AA.inversed();
    auto t7 = std::chrono::high_resolution_clock::now();
    MatrixXd invA3=AA.inversed_rob();
    auto t8 = std::chrono::high_resolution_clock::now();
    std::cout <<"Inversa: "<< std::chrono::duration_cast<std::chrono::microseconds>(t6-t5).count()<<"\n";
    std::cout <<"Inversa LU: "<< std::chrono::duration_cast<std::chrono::microseconds>(t7-t6).count()<<"\n";
    std::cout <<"Inversa LUP: "<< std::chrono::duration_cast<std::chrono::microseconds>(t8-t7).count()<<"\n";
    std::cout <<"\nInversa Gauss:\n";
    mprint(invA1);
    std::cout <<"\nInversa LU:\n";
    mprint(invA2);
    std::cout <<"\nInversa LUP:\n";
    mprint(invA3);
    // TEST FATTORIZZAZIONE
    std::cout << "\nTest Fattorizzazione";
    MatrixXd LL(AA);
    MatrixXd UU(AA);
    MatrixXs PP(4,1);
    auto t9 = std::chrono::high_resolution_clock::now();
    LU_Crout(AA, LL, UU);
    auto t10 = std::chrono::high_resolution_clock::now();
    LU_Cormen(AA, LL, UU);
    auto t11 = std::chrono::high_resolution_clock::now();
    LUP_Cormen(AA, LL, UU, PP);
    auto t12 = std::chrono::high_resolution_clock::now();
    std::cout <<"\nFatt. Crout LU: "<< std::chrono::duration_cast<std::chrono::microseconds>(t10-t9).count()<<"\n";
    std::cout <<"Fatt. Cormen LU: "<< std::chrono::duration_cast<std::chrono::microseconds>(t11-t10).count()<<"\n";
    std::cout <<"Fatt. Cormen LUP: "<< std::chrono::duration_cast<std::chrono::microseconds>(t12-t11).count()<<"\n";
    printf("\nDeterminante di AA: %f\n",AA.det());
    // TEST QUAD
    std::cout << "\nTest Quad_mult"<<"\n";
    auto t13 = std::chrono::high_resolution_clock::now();
    MatrixXd RES1=QuadProd(AA, QuadProd((~AA), !BB));
    auto t14 = std::chrono::high_resolution_clock::now();
    MatrixXd RES2=AA*~AA*(!BB)*AA*(~AA);
    auto t15 = std::chrono::high_resolution_clock::now();
    std::cout <<"Funzione Nuova: "<< std::chrono::duration_cast<std::chrono::microseconds>(t14-t13).count()<<"\n";
    std::cout <<"Metodo Standard: "<< std::chrono::duration_cast<std::chrono::microseconds>(t15-t14).count()<<"\n";
    std::cout << "\nMetodo mio:"<<"\n";
    mprint(RES1);
    std::cout << "\nMetodo standard:"<<"\n";
    mprint(RES2);
    return 0;
}