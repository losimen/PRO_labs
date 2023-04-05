#include <vector>
#include <random>
#include <iostream>


typedef std::vector<std::vector<double>> Matrix;
typedef std::vector<double> Vector;


namespace MatrixUtils
{
    Matrix Transpose(const Matrix &matrix)
    {
        Matrix result;
        result.reserve(matrix.size());

        for (unsigned i = 0; i < matrix.size(); ++i)
        {
            Vector row;
            row.reserve(matrix.size());

            for (const auto & j : matrix)
            {
                row.emplace_back(j[i]);
            }

            result.emplace_back(row);
        }

        return result;
    }

    Matrix Multiply(const Matrix &matrix1, const Matrix &matrix2)
    {
        Matrix result;
        result.reserve(matrix1.size());

        for (const auto & i : matrix1)
        {
            Vector row;
            row.reserve(i.size());

            for (unsigned j = 0; j < i.size(); ++j)
            {
                double sum = 0;
                for (unsigned k = 0; k < i.size(); ++k)
                {
                    sum += i[k] * matrix2[k][j];
                }

                row.emplace_back(sum);
            }

            result.emplace_back(row);
        }

        return result;
    }

    Matrix Sub(const Matrix &matrix1, const Matrix &matrix2)
    {
        Matrix result;
        result.reserve(matrix1.size());

        for (unsigned i = 0; i < matrix1.size(); ++i)
        {
            Vector row;
            row.reserve(matrix1[i].size());

            for (unsigned j = 0; j < matrix1[i].size(); ++j)
            {
                row.emplace_back(matrix1[i][j] - matrix2[i][j]);
            }

            result.emplace_back(row);
        }

        return result;
    }
}


void printMatrix(const Matrix &matrix)
{
    for (const auto &row: matrix)
    {
        for (const auto &elem: row)
        {
            std::cout << elem << " ";
        }
        std::cout << std::endl;
    }
}


void randGen(Matrix &A, int N, int start = 1, int end = 100)
{
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(start, end);

    A.reserve(N);
    for (unsigned i = 0; i < N; ++i)
    {
        Vector row;
        row.reserve(N);

        for (unsigned j = 0; j < N; ++j)
        {
            row.emplace_back((int)dist(mt));
        }

        A.emplace_back(row);
    }
}


class Phase1
{
public:
    const int N;
    Vector vectorColumnBi {};
    Matrix A {};

    explicit Phase1(int n)
            : N(n)
    {}

    Vector countY1()
    {
        Vector result;
        result.reserve(vectorColumnBi.size());

        for (const auto &row: A)
        {
            double sum = 0;
            for (unsigned i = 0; i < row.size(); ++i)
            {
                sum += row[i] * vectorColumnBi[i];
            }
            result.emplace_back(sum);
        }

        return result;
    }

    void genBi(int num)
    {
        vectorColumnBi.reserve(N);

        for (unsigned i = 0; i < N; ++i)
        {
            vectorColumnBi.emplace_back(int(num * (i + 1)));
        }
    }

    void genA()
    {
        randGen(A, N);
    }
};


class Phase2
{
public:
    const int N;
    Matrix A1;
    Vector vector_column_bi;
    Vector vector_column_ci;

    explicit Phase2(int n)
            : N(n)
    {}

    void genBi(int num)
    {
        vector_column_bi.reserve(N);

        for (unsigned i = 0; i < N; ++i)
        {
            vector_column_bi.emplace_back(int(num * (i + 1)));
        }
    }

    void genCi(int num)
    {
        vector_column_ci.reserve(N);

        for (unsigned i = 0; i < N; ++i)
        {
            vector_column_ci.emplace_back(int(num * (i + 1)));
        }
    }

    void genA1()
    {
        randGen(A1, N);
    }

    Vector countY2()
    {
        auto subTwoVectors = [](const Vector &v1, const Vector &v2)
        {
            std::vector<int> result;
            result.reserve(v1.size());

            for (unsigned i = 0; i < v1.size(); ++i)
            {
                result.emplace_back(v1[i] - v2[i]);
            }

            return result;
        };

        Vector result;
        result.reserve(vector_column_ci.size());

        for (const auto &row: A1)
        {
            const auto sub = subTwoVectors(row, vector_column_bi);

            double sum = 0;
            for (unsigned i = 0; i < sub.size(); ++i)
            {
                sum += sub[i] * vector_column_ci[i];
            }

            result.emplace_back(sum);
        }

        return result;
    }
};


class Phase3
{
public:
    const int N;
    Matrix A2;
    Matrix B2;
    Matrix C2;

    explicit Phase3(int n)
            : N(n)
    {}

    void genA2(int start = 1, int end = 100)
    {
        randGen(A2, N, start, end);
    }

    void genB2(int start = 1, int end = 100)
    {
        randGen(B2, N, start, end);
    }

    void genC2()
    {
        C2.reserve(N);
        for (int i = 0; i < N; ++i)
        {
            C2.reserve(N);
            C2.emplace_back();
            for (int j = 0; j < N; ++j)
            {
                C2[i].emplace_back(1 / (double)(i + 2*j+1));
            }
        }

        std::cout << "C2: " << C2.size() << std::endl;
    }

    Matrix countY3() const
    {
        // A2 * (B2 - C2)
        const auto sub = MatrixUtils::Sub(B2, C2);
        return MatrixUtils::Multiply(A2, sub);
    }
};


int main()
{
    const int N = 4;

    Phase1 p1(N);
    p1.genA();
    p1.genBi(9);
    Vector y1 = p1.countY1();

    Matrix y1_transposed = MatrixUtils::Transpose({y1});
    printMatrix(y1_transposed);

    Phase2 p2(N);
    p2.genA1();
    p2.genBi(9);
    p2.genCi(7);
    Vector y2 = p2.countY2();

    Phase3 p3(N);
    p3.genA2(1, 15);
    p3.genB2(1, 15);
    p3.genC2();
    Matrix Y3 = p3.countY3();
    Matrix Y3_2 = MatrixUtils::Multiply(Y3, Y3);
    Matrix Y3_3 = MatrixUtils::Multiply(Y3_2, Y3);

    return 0;
}
