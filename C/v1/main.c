//
//  main.c
//  asdasd
//
//  Created by Andrea Vivani on 14/08/16.
//  Copyright © 2016 Andrea Vivani. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Matrix.h"
#include "NumMethods.h"

typedef union{
    uint16_t ALL;
    struct {
        uint8_t W1 : 1;
        uint8_t W2 : 1;
        uint8_t W3 : 1;
        uint8_t W4 : 1;
        uint8_t W5 : 1;
        uint8_t W6 : 1;
        uint8_t W7 : 1;
        uint8_t W8 : 1;
    };
}warnings_t;

struct{
    uint8_t bb [3];
    uint8_t aa[4];
}test;

void mprint(Matrix m){
    printf("%dx%d\n", m->rows, m->cols);
    for(int i=0; i<m->rows; i++) {
        for (int j = 0; j<m->cols; j++)
            printf("%4.6f\t",(float) m->data[(i) * m->cols + (j)]);
        printf("\n");
    }
}


int main(int argc, const char * argv[]) {
    // TEST GAUSS
    printf("Test Gauss\n");
    /*float DataA []={-0.009345794,	0.00952381,	1,
     -1,	-0.00952381,	0.044247788,
     1,	0.028571429,	0.008849558,
     -0.009345794,	1,	-0.115044248,
     0.028037383,	-1,	0.008849558,
     -0.775700935,	-0.638095238,	0.008849558,
     -0.570093458,	0.80952381,	-0.026548673,
     0.644859813,	0.771428571,	-0.061946903,
     0.775700935,	-0.619047619,	-0.115044248,
     0.981308411,	0.00952381,	-0.256637168};*/
    
    /*float DataA []={0.6568, -0.6316, -9.7223,
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
    
    /*float DataA []={-0.054563, 0.140088, -10.569158,
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
    
    float DataA []={0.207943, 0.176336, -10.472851,
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
    Matrix Data = newMatrixData(20, 3, DataA);
    float _data0[]={0,0,0,1,0,0,1,0,1};
    Matrix X0 = newMatrixData(9, 1, _data0);
    float _data01[]={0,0,0,1,1,1};
    Matrix X01 = newMatrixData(6, 1, _data01);
    Matrix P = GaussNewton_Sens_Cal_9(Data, 9.81, X0, 600, 1e-9);
    Matrix P2 = GaussNewton_Sens_Cal_6(Data, 9.81, X01, 600, 1e-9);
    printf("\nData\n");
    mprint(Data);
    printf("\n9 Parametri:\n");
    mprint(P);
    printf("\n6 Parametri:\n");
    mprint(P2);
    matFree(Data);
    matFree(X0);
    matFree(X01);
    matFree(P);
    matFree(P2);
    // TEST LINSOLVE
    printf("\nTest LinSolve\n");
    float _dataA[]={0.5432,    0.3171,    0.3816,    0.4898,
        0.0462,    0.4358,    0.6651,    0.4456,
        0.8235,    0.1324,    0.7952,    0.6463,
        0.6948,    0.9745,    0.1869,    0.4456};
    float _dataB[]={0.7547,    0.1626,    0.3404,    0.2551,
        0.2760,    0.1190,    0.5853,    0.5060,
        0.6797,    0.4984,    0.2238,    0.6991,
        0.6551,    0.9597,    0.7513,   0.8909};
    Matrix AA = newMatrixData(4, 4, _dataA);
    Matrix BB = newMatrixData(4, 4, _dataB);
    Matrix tmp1 = newMatrix(1, 1);
    matInversed(AA, tmp1);
    Matrix Soli = newMatrix(1, 1);
    matMult(tmp1, BB, Soli);
    printf("\nSoluz. Inversa:\n");
    mprint(Soli);
    Matrix Solg = newMatrix(1,1);
    LinSolveGauss(AA, BB, Solg);
    printf("\nSoluz. Gauss:\n");
    mprint(Solg);
    Matrix Soll = newMatrix(1,1);
    LinSolveLU(AA, BB, Soll);
    printf("\nSoluz. LU:\n");
    mprint(Soll);
    Matrix Soll2 = newMatrix(1,1);
    LinSolveLUP(AA, BB, Soll2);
    printf("\nSoluz. LUP:\n");
    mprint(Soll2);
    matFree(tmp1);
    matFree(Soli);
    matFree(Solg);
    matFree(Soll);
    matFree(Soll2);
    // TEST INVERSA
    printf("\nTest Inversa\n");
    Matrix invA1 = newMatrix(1, 1);
    Matrix invA2 = newMatrix(1, 1);
    matInversed(AA, invA1);
    matInversed_rob(AA, invA2);
    printf("\nInversa LU:\n");
    mprint(invA1);
    printf("\nInversa LUP:\n");
    mprint(invA2);
    matFree(invA1);
    matFree(invA2);
    // TEST FATTORIZZAZIONE
    printf("\nTest Fattorizzazione");
    Matrix LL = copyMatrix(AA);
    Matrix UU = copyMatrix(AA);
    Matrix PP = newMatrix(4,1);
    LU_Crout(AA, LL, UU);
    LU_Cormen(AA, LL, UU);
    LUP_Cormen(AA, LL, UU, PP);
    printf("\nDeterminante di AA: %f\n", matDet(AA));
    matFree(LL);
    matFree(UU);
    matFree(PP);
    // TEST QUAD
    printf("\nTest Quad_mult\n");
    tmp1 = newMatrix(1, 1);
    Matrix tmp2 = newMatrix(1, 1);
    Matrix tmp3 = newMatrix(1, 1);
    Matrix tmp4 = newMatrix(1, 1);
    Matrix RES1 = newMatrix(1, 1);
    Matrix RES2 = newMatrix(1, 1);
    matInversed(BB, tmp1);
    matTrans(AA, tmp2);
    QuadProd(tmp2, tmp1, tmp3); //tmp3 = ~AA * !BB * AA
    QuadProd(AA, tmp3, RES1); // RES1 = AA * ~AA * !BB * AA * ~AA
    //
    matMult(AA, tmp2, tmp3); // tmp3 = AA * ~AA
    matMult(tmp3, tmp1, tmp4); // tmp4 = AA * ~AA * !BB
    matMult(tmp4, AA, tmp3); // tmp3 = AA * ~AA * !BB * AA
    matMult(tmp3, tmp2, RES2); // RES2 = AA * ~AA * !BB * AA * ~AA
    printf("\nMetodo mio:\n");
    mprint(RES1);
    printf("\nMetodo standard:\n");
    mprint(RES2);
    matFree(AA);
    matFree(BB);
    matFree(tmp1);
    matFree(tmp2);
    matFree(tmp3);
    matFree(tmp4);
    matFree(RES1);
    matFree(RES2);
    
    warnings_t warnings;
    warnings.W8 = 0;
    warnings.W5 = 0;
    warnings.W2 = 1;
    warnings.W1 = 0;
    printf("\r\n%x\r\n", warnings.ALL);
    return 0;
}
