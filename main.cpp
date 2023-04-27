#include <random>
#include <iostream>
#include <iomanip>
#include <vector>

unsigned SIZE_N = 10;

typedef long double CalVar;
typedef std::vector<CalVar> TVector;
typedef std::vector<TVector> TMatrix;
typedef std::vector<TMatrix> TParallelMatrix;


CalVar randGenNumber(int start = 1, int end = 100)
{
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<CalVar> dist(start, end);

    return dist(mt);
}


void genMatrixB(TMatrix &B)
{
    for(int i = 0; i < SIZE_N; i++) {
        for(int j = 0; j < SIZE_N; j++) {
            if(i < SIZE_N / 2)
            {
                if(j < i || j >= SIZE_N - i)
                    B[i][j] = 0;
                else
                    B[i][j] = randGenNumber(1, 100);
            }
            else {
                if(j < SIZE_N - i - 1 || j > i)
                    B[i][j] = 0;
                else
                    B[i][j] = randGenNumber(1, 100);
            }
        }
    }
}


void genMatrixA(TMatrix &A)
{
    for(int i = 0; i < SIZE_N; i++) {
        for(int j = 0; j < SIZE_N; j++) {
            if(j < SIZE_N - i)
                A[i][j] = i + 1;
            else
                A[i][j] = 0;
        }
    }
}


void initMatrix(TMatrix &matrix)
{
    matrix.reserve(SIZE_N);

    for (int i = 0; i < SIZE_N; i++)
    {
        matrix[i].reserve(SIZE_N);
    }
}


void initParallelMatrix(TParallelMatrix &matrix)
{
    matrix.reserve(SIZE_N);

    for (int i = 0; i < SIZE_N; i++)
    {
        matrix[i].reserve(SIZE_N);

        for (int j = 0; j < SIZE_N; j++)
        {
            matrix[i][j].reserve(SIZE_N);
        }
    }
}


void printMatrix(TMatrix &matrix)
{
    for(int i = 0; i < SIZE_N; i++) {
        for(int j = 0; j < SIZE_N; j++) {
            std::cout << std::setw(5) << matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
}


void printParallelMatrix(TParallelMatrix &matrix)
{
    for(int i = 0; i < SIZE_N; i++) {
        for(int j = 0; j < SIZE_N; j++) {
            std::cout << std::setw(5) << matrix[i][j][SIZE_N-1] << " ";
        }
        std::cout << std::endl;
    }
}


unsigned mulStandard(TParallelMatrix &result, TMatrix &matrix1, TMatrix &matrix2)
{
    unsigned numberOfOperations = 0;

    for (size_t i = 0; i < SIZE_N; i++) {
        for (size_t j = 0; j < SIZE_N; j++) {
            result[i][j][0] = 0;
            for (size_t k = 0; k < SIZE_N; k++) {
                result[i][j][k+1] = result[i][j][k] + matrix1[i][k] * matrix2[k][j];
                numberOfOperations += 2;
            }
        }
    }

    return numberOfOperations;
}


unsigned mulOptimisedRecursive(TParallelMatrix &result, TMatrix &matrix1, TMatrix &matrix2)
{
    static unsigned numberOfOperations = 0;

    static unsigned i = 0;
    static unsigned j = 0;
    static unsigned k = 0;

    if (i < SIZE_N)
    {
        if (j < SIZE_N)
        {
            result[i][j][0] = 0;
            if (k < SIZE_N)
            {
                if ( k >= SIZE_N - 1 )
                {
                    result[i][j][k+1] = result[i][j][k];
                    k++;
                    mulOptimisedRecursive(result, matrix1, matrix2);
                    return 0;
                }

                if ( k < SIZE_N / 2)
                {
                    if( j < k || j >= SIZE_N - k )
                    {
                        result[i][j][k+1] = result[i][j][k];
                        k++;
                        mulOptimisedRecursive(result, matrix1, matrix2);
                        return 0;
                    }
                }
                else
                {
                    if ( j < SIZE_N - k - 1 || j > k)
                    {
                        result[i][j][k+1] = result[i][j][k];
                        k++;
                        mulOptimisedRecursive(result, matrix1, matrix2);
                        return 0;
                    }
                }

                result[i][j][k+1] = result[i][j][k] + matrix1[i][k] * matrix2[k][j];
                numberOfOperations += 2;
                k++;
                mulOptimisedRecursive(result, matrix1, matrix2);
            }
            k = 0;
            j++;
            mulOptimisedRecursive(result, matrix1, matrix2);
        }
        j = 0;
        i++;
        mulOptimisedRecursive(result, matrix1, matrix2);
    }

    return numberOfOperations;
}



int main()
{
    TMatrix B;
    initMatrix(B);
    genMatrixB(B);

    TMatrix A;
    initMatrix(A);
    genMatrixA(A);

    TParallelMatrix resultStandard;
    initParallelMatrix(resultStandard);
    unsigned numberOfOperationsStandard = mulStandard(resultStandard,A, B);

    TParallelMatrix resultParallelRecursive;
    initParallelMatrix(resultParallelRecursive);
    unsigned numberOfOperationsRecursive = mulOptimisedRecursive(resultParallelRecursive, A, B);

    std::cout << std::endl << "Result standard:" << std::endl;
    printParallelMatrix(resultStandard);
    std::cout << "Number of operations standard: " << numberOfOperationsStandard << std::endl;


    std::cout << std::endl << "Result parallel: " << std::endl;
    printParallelMatrix(resultParallelRecursive);
    std::cout << "Number of operations recursive: " << numberOfOperationsRecursive << std::endl;

    return 0;
}