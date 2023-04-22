#include <iostream>
#include <array>

constexpr int SIZE_N = 6;
typedef std::array<std::array<long double, SIZE_N>, SIZE_N> Matrix;


Matrix genMatrixA()
{
    // 111…111
    //222…220
    //333…300
    //…
    //n000…00

    Matrix A;
    for (int i = 0; i < SIZE_N; ++i)
    {
        for (int j = 0; j < SIZE_N; ++j)
        {
            A[i][j] = i * 10 + j;
        }
    }

    return A;
}


int main()
{


    return 0;
}
