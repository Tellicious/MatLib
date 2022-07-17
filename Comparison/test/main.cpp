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
#include "MatrixF.h"
#include "NumMethodsF.h"

template <typename T> void mprint(MatrixX<T>& r) {
    std::cout << "\n" << (int) r.rows() << "x" << (int) r.columns() << "\n";
    for(int i=0; i<r.rows(); i++) {
        for (int j=0; j<r.columns(); j++)
            printf("%4.6f\t",(double) r(i,j));
        std::cout << "\n";
    }
}

void mprint2(MatrixF m){
    std::cout << "\n" << (int) m.rows << "x" << (int) m.cols << "\n";
    for(int i=0; i<m.rows; i++) {
        for (int j=0; j<m.cols; j++)
            printf("%4.6f\t",(float) matGet(m,i,j));
        std::cout << "\n";
    }

}

float asdasd(float a, float b=3){
    return a*b;
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

class Test{
public:
    Test():M1(3,1){
        instanced++;
        instance=instanced;};
    static int instanced;
    int instance;
    uint8_t setta(){
        M1(0,0)=33;
        M1(1,0)=12;
        M1(2,0)=13;
        return 1;
    }
    MatrixXd M1;
           
};
int Test::instanced=0;

float constrain(float x,float inf,float sup){
if (x<inf)
return inf;
else if (x>sup)
return sup;
else
return x;
}

MatrixXf motor_mix_4_lin(MatrixXf U, float omega_0, float omega_min, float omega_max){
    
    omega_0 = fmaxf(omega_0, omega_min);
    
    U(0,0) = 4 * constrain(0.25 * U(0,0), (omega_min - omega_0), (omega_max - omega_0));
    float omega_th = omega_0 + 0.25 * U(0,0);
    
    float Du = 2 * fminf(fabsf(omega_max - omega_th), fabsf(omega_min - omega_th));
    U(1,0) = constrain(U(1,0), - Du, Du);
    U(2,0) = constrain(U(2,0), - Du, Du);
    
    float U_m = fmaxf(fabsf(U(1,0)), fabsf(U(2,0)));
    Du = 4 * fminf(fabs(omega_max - omega_th - 0.5 * U_m), fabs(omega_min - omega_th + 0.5 * U_m));
    U(3,0) = constrain(U(3,0), -Du , Du);
    
    MatrixXf M(4,1);
    M(0,0) = omega_th + 0.5 * U(2,0) - 0.25 * U(3,0);
    M(1,0) = omega_th - 0.5 * U(1,0) + 0.25 * U(3,0);
    M(2,0) = omega_th - 0.5 * U(2,0) - 0.25 * U(3,0);
    M(3,0) = omega_th + 0.5 * U(1,0) + 0.25 * U(3,0);
    //fprintf('m1=%f, m2=%f, m3=%f, m4=%f\n',m1,m2,m3,m4);
    //fprintf('U1=%f, U2=%f, U3=%f, U4=%f\n',m1+m2+m3+m4-4*om0,m4-m2,m1-m3,m2+m4-m1-m3);
    return M;
};

int main(int argc, const char * argv[]) {
    // TEST GAUSS
    std::cout<<"Test Gauss\n";
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
    
    /*Data={0.6568, -0.6316, -9.7223,
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
        -0.1903, -0.1735, -9.8529};*/
    
    /*Data={-0.054563, 0.140088, -10.569158,
        1.467040, 1.883719, 8.909904,
        1.339795, -1.834797, 8.971394,
        -2.660742, -1.043005, 8.801044,
        -2.180823, 1.632330, 8.898216,
        0.241814, 2.300491, -10.292187,
        -1.895915, -0.029317, -10.336355,
        2.104321, -2.263695, -10.034464,
        2.003852, 2.056409, -10.095016,
        -1.788120, 2.208810, -10.153329,
        0.191610, 9.830813, -0.349682,
        9.751721, -0.103503, -0.964587,
        -0.115117, -9.703846, -0.427194,
        -9.835476, 0.154282, -0.969580,
        -6.964331, 6.979779, -0.881884,
        6.967795, 6.900365, -0.751046,
        6.876878, -6.856656, -0.728948,
        -6.935593, -6.858181, -1.142173, 
        -6.874837, 7.064003, -0.966496, 
        0.764942, -0.097632, 9.206914};*/
    
    double Dataaa[]={0.207943, 0.176336, -10.472851,
        -0.032664, -9.691158, -0.525727,
        9.739063, 0.100049, -0.033633,
        -0.014726, 9.848834, -0.441255,
        -9.848631, 0.111619, -0.447592,
        -0.817354, 1.245916, -10.409218,
        -0.954637, -0.850694, -10.362207,
        1.181917, -0.940953, -10.381795,
        1.234348, 1.037082, -10.323170,
        0.493486, 2.144262, 9.028279,
        -1.948339, 1.323566, 9.077987,
        -1.715729, -1.417367, 9.021542,
        1.402404, -1.214459, 9.140464,
        6.862244, 7.007604, -0.370883,
        7.013921, -6.699038, -0.329180,
        -6.915473, -6.891760, -0.527419,
        -7.064906, 6.897670, -0.619903,
        3.656467, 3.660516, -9.034408, 
        3.640905, -3.403118, -9.039276, 
        -3.512549, -3.599191, -9.047305};
    float Dataaa2[]={0.207943, 0.176336, -10.472851,
        -0.032664, -9.691158, -0.525727,
        9.739063, 0.100049, -0.033633,
        -0.014726, 9.848834, -0.441255,
        -9.848631, 0.111619, -0.447592,
        -0.817354, 1.245916, -10.409218,
        -0.954637, -0.850694, -10.362207,
        1.181917, -0.940953, -10.381795,
        1.234348, 1.037082, -10.323170,
        0.493486, 2.144262, 9.028279,
        -1.948339, 1.323566, 9.077987,
        -1.715729, -1.417367, 9.021542,
        1.402404, -1.214459, 9.140464,
        6.862244, 7.007604, -0.370883,
        7.013921, -6.699038, -0.329180,
        -6.915473, -6.891760, -0.527419,
        -7.064906, 6.897670, -0.619903,
        3.656467, 3.660516, -9.034408,
        3.640905, -3.403118, -9.039276,
        -3.512549, -3.599191, -9.047305};
    MatrixXd Data(20,3,Dataaa);
    MatrixF Datav2 = newMatrixData(20, 3, Dataaa2);
    double _data0[]={0,0,0,1,0,0,1,0,1};
    float _data0v2[]={0,0,0,1,0,0,1,0,1};
    MatrixXd X0(9,1,_data0);
    MatrixF X0v2 = newMatrixData(9, 1, _data0v2);
    double _data01[]={0,0,0,1,1,1};
    float _data01v2[]={0,0,0,1,1,1};
    MatrixXd X01(6,1,_data01);
    MatrixF X01v2 = newMatrixData(6, 1, _data01v2);
    MatrixXd P=GaussNewton_Sens_Cal_9(Data,9.81,X0, 600,1e-9);
    MatrixXd P2=GaussNewton_Sens_Cal_6(Data, 9.81,X01, 600, 1e-9);
    MatrixF Pv2=GaussNewton_Sens_Cal_9F(Datav2, 9.81, X0v2, 600, 1e-9);
    MatrixF P2v2=GaussNewton_Sens_Cal_6F(Datav2, 9.81, X01v2, 600, 1e-9);
    std::cout<<"\nData\n";
    mprint(Data);
    std::cout<<"\n9 Parametri:\n";
    mprint(P);
    mprint2(Pv2);
    std::cout<<"\n6 Parametri:\n";
    mprint(P2);
    mprint2(P2v2);
    // TEST LINSOLVE
    std::cout<<"\nTest LinSolve\n";
    double _dataA[]={0.5432,    0.3171,    0.3816,    0.4898,
        0.0462,    0.4358,    0.6651,    0.4456,
        0.8235,    0.1324,    0.7952,    0.6463,
        0.6948,    0.9745,    0.1869,    0.4456};
    float _dataA2[]={0.5432,    0.3171,    0.3816,    0.4898,
        0.0462,    0.4358,    0.6651,    0.4456,
        0.8235,    0.1324,    0.7952,    0.6463,
        0.6948,    0.9745,    0.1869,    0.4456};
    double _dataB[]={0.7547,    0.1626,    0.3404,    0.2551,
        0.2760,    0.1190,    0.5853,    0.5060,
        0.6797,    0.4984,    0.2238,    0.6991,
        0.6551,    0.9597,    0.7513,   0.8909};
    float _dataB2[]={0.7547,    0.1626,    0.3404,    0.2551,
        0.2760,    0.1190,    0.5853,    0.5060,
        0.6797,    0.4984,    0.2238,    0.6991,
        0.6551,    0.9597,    0.7513,   0.8909};
    MatrixXd AA(4,4,_dataA);
    MatrixXd BB(4,4,_dataB);
    MatrixF AA2 = newMatrixData(4, 4, _dataA2);
    MatrixF BB2 = newMatrixData(4, 4, _dataB2);
    auto t0 = std::chrono::high_resolution_clock::now();
    MatrixXd Soli=!AA*BB;
    //MatrixF Soli2 = matMult(matInversed(AA2), BB2);
    auto t1 = std::chrono::high_resolution_clock::now();
    MatrixXd Solg=LinSolveGauss(AA, BB);
    MatrixF Solg2 = LinSolveGaussF(AA2, BB2);
    auto t2 = std::chrono::high_resolution_clock::now();
    MatrixXd Soll=LinSolveLU(AA, BB);
    MatrixF Sollv2=LinSolveLUF(AA2, BB2);
    auto t3 = std::chrono::high_resolution_clock::now();
    MatrixXd Soll2=LinSolveLUP(AA, BB);
    MatrixF Soll2v2=LinSolveLUPF(AA2, BB2);
    auto t4 = std::chrono::high_resolution_clock::now();
    std::cout <<"Inversa: "<< std::chrono::duration_cast<std::chrono::microseconds>(t1-t0).count()<<"\n";
    std::cout <<"Gauss: "<< std::chrono::duration_cast<std::chrono::microseconds>(t2-t1).count()<<"\n";
    std::cout <<"LU: "<< std::chrono::duration_cast<std::chrono::microseconds>(t3-t2).count()<<"\n";
    std::cout <<"LUP: "<< std::chrono::duration_cast<std::chrono::microseconds>(t4-t3).count()<<"\n";
    std::cout <<"\nSoluz. Inversa:\n";
    mprint(Soli);
    //mprint2(Soli2);
    std::cout <<"\nSoluz. Gauss:\n";
    mprint(Solg);
    mprint2(Solg2);
    std::cout <<"\nSoluz. LU:\n";
    mprint(Soll);
    mprint2(Sollv2);
    std::cout <<"\nSoluz. LUP:\n";
    mprint(Soll2);
    mprint2(Soll2v2);
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
    MatrixXd RES2(1,1);
    RES2=AA*(~AA)*(!BB)*AA*(~AA);
    auto t15 = std::chrono::high_resolution_clock::now();
    std::cout <<"Funzione Nuova: "<< std::chrono::duration_cast<std::chrono::microseconds>(t14-t13).count()<<"\n";
    std::cout <<"Metodo Standard: "<< std::chrono::duration_cast<std::chrono::microseconds>(t15-t14).count()<<"\n";
    std::cout << "\nMetodo mio:"<<"\n";
    mprint(RES1);
    std::cout << "\nMetodo standard:"<<"\n";
    mprint(RES2);
    Test aa;
    Test bb;
    Test cc;
    aa.setta();
    mprint(aa.M1);
    std::cout<<Test::instanced<<"\n"<<aa.instance<<"\n"<<bb.instance<<"\n"<<cc.instance<<"\n";
    std::cout<<asdasd(4,5)<<","<<asdasd(2)<<"\n";
    MatrixXf U(4,1);
    U(0,0)=50;
    U(1,0)=10;
    U(2,0)=-2;
    U(3,0)=-20;
    MatrixXf M = motor_mix_4_lin(U,50,120,200);
    mprint(M);
    
    std::cout<<"\n\n\n\n\n";
    float dd1[9]={3,4,1,9,5,3,6,8,4};
    MatrixF asd1 = newMatrixData(3, 3, dd1);
    mprint2(asd1);
    asd1 = matSubScalar(asd1, 4);
    mprint2(asd1);
    mprint2(matAdd(asd1, asd1));
    mprint2(matMult(asd1, asd1));
    
    
    uint64_t address = 0x0102030405;
    uint8_t add[5];
    memcpy(add, &address, sizeof(add));
    std::cout<<(int) add[0]<<","<<(int) add[1]<<","<<(int) add[2]<<","<<(int) add[3]<<","<<(int) add[4]<<"\n";
    return 0;
}