# MatLib
### C++ lightweight Matrix Library with several numeric methods included (system solvers, least squares fittings, etc.)

- Include with `#include "MatLib.h"`

- Available types:
    ```cpp
    typedef MatrixX<int8_t>      MatrixXs;
    typedef MatrixX<uint8_t>     MatrixXus;
    typedef MatrixX<int16_t>     MatrixXi;
    typedef MatrixX<uint16_t>    MatrixXui;
    typedef MatrixX<int32_t>     MatrixXl;
    typedef MatrixX<uint32_t>    MatrixXul;
    typedef MatrixX<float>       MatrixXf;
    typedef MatrixX<double>      MatrixXd;
    ```

- Example of constructor with data:
    ```cpp
    double _data0[] = {0, 0, 0, 1, 0, 0, 1, 0, 1};
    MatrixXd X0(9, 1, _data0);
    ```

- Example of constructor without data:
    ```cpp
    MatrixXf AA(3, 3);
    MatrixXf BB(4, 4);
    ```

- Example of usage:
    ```cpp
    MatrixXf RES2 = AA * (~AA) * (!BB) * AA * (~AA);
    ```

- Description of numeric methods included available in `Num_methods.h`


