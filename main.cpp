#include <random>
#include <iostream>
#include <array>
#include <iomanip>

constexpr unsigned SIZE_N = 4;

template <unsigned N>
class Matrix
{
private:
    typedef std::array<std::array<double, N>, N> TMatrix;
    typedef std::array<double, N> TVector;
    TMatrix m_matrix;

public:
    Matrix() = delete;
    Matrix(unsigned rows, unsigned cols): rows(rows), cols(cols)
    {
        m_matrix.fill(TVector());
    }

    TVector& operator[](unsigned index) {
        return m_matrix[index];
    }

    unsigned rows;
    unsigned cols;

    unsigned size() const
    {
        return m_matrix.size();
    }

    Matrix transpose()
    {
        Matrix<SIZE_N> result(rows, cols);

        for (unsigned i = 0; i < m_matrix.size(); ++i)
            for (unsigned j = 0; j < m_matrix[i].size(); ++j)
                result[i][j] = m_matrix[j][i];

        return result;
    }

    Matrix mul(Matrix &oth)
    {
        if (m_matrix.size() != oth[0].size())
            throw std::invalid_argument("m_matrix sizes are not equal");

        Matrix<SIZE_N> result(m_matrix.size(), oth[0].size());

        for (unsigned i = 0; i < m_matrix.size(); ++i)
        {
            for (unsigned j = 0; j < m_matrix[i].size(); ++j)
            {
                double sum = 0;
                for (unsigned k = 0; k < m_matrix[i].size(); ++k)
                {
                    sum += m_matrix[i][k] * oth[k][j];
                }

                result[i][j] = sum;
            }
        }

        return result;
    }

    Matrix sub(Matrix &oth)
    {
        if (m_matrix.size() != oth.size())
            throw std::invalid_argument("m_matrix sizes are not equal");

        Matrix<SIZE_N> result(m_matrix.size(), oth[0].size());

        for (unsigned i = 0; i < m_matrix.size(); ++i)
        {
            for (unsigned j = 0; j < m_matrix[i].size(); ++j)
            {
                result[i][j] = m_matrix[i][j] - oth[i][j];
            }
        }

        return result;
    }

    Matrix add(Matrix &oth)
    {
        if (m_matrix.size() != oth.size())
            throw std::invalid_argument("m_matrix sizes are not equal");

        Matrix<SIZE_N> result(m_matrix.size(), oth[0].size());

        for (unsigned i = 0; i < m_matrix.size(); ++i)
        {
            for (unsigned j = 0; j < m_matrix[i].size(); ++j)
            {
                result[i][j] = m_matrix[i][j] + oth[i][j];
            }
        }

        return result;
    }

    void printMatrix()
    {
        for (unsigned i = 0; i < m_matrix.size(); ++i)
        {
            for (unsigned j = 0; j < m_matrix[i].size(); ++j)
            {
                std::cout << std::setw(10) << std::setprecision(5) << m_matrix[i][j] << " ";
            }
            std::cout << std::endl;
        }
    }
};

void randGenMatrix(Matrix<SIZE_N> &A, int N, int start = 1, int end = 100)
{
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(start, end);

    for (unsigned i = 0; i < N; ++i)
    {
        for (unsigned j = 0; j < N; ++j)
        {
            A[i][j] = dist(mt);
        }
    }
}

void randGenColumnVector(Matrix<SIZE_N> &A, int N, int start = 1, int end = 100)
{
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(start, end);

    for (unsigned i = 0; i < N; ++i)
    {
        A[i][0] = dist(mt);
    }
}


class Phase1
{
public:
    const int N;
    Matrix<SIZE_N> vectorColumnBi = Matrix<SIZE_N>(N, 1);
    Matrix<SIZE_N> A = Matrix<SIZE_N>(N, N);

    explicit Phase1(int n)
            : N(n)
    {
        genBi(9);
        genA();
    }

    Matrix<SIZE_N> countY1()
    {
        Matrix<SIZE_N> result(N, 1);

        for (unsigned i = 0; i < N; ++i)
        {
            double sum = 0;
            for (unsigned j = 0; j < N; ++j)
            {
                sum += A[i][j] * vectorColumnBi[i][0];
            }
            result[i][0] = sum;
        }

        return result;
    }

    void genBi(int num)
    {
        for (unsigned i = 0; i < N; ++i)
        {
            vectorColumnBi[i][0] = int(num * (i + 1));
        }
    }

    void genA(int start = 0, int end = 100)
    {
        randGenMatrix(A, N, start, end);
    }
};


class Phase2
{
public:
    const int N;
    Matrix<SIZE_N> A1 = Matrix<SIZE_N>(N, N);
    Matrix<SIZE_N> vectorColumnBi = Matrix<SIZE_N>(N, 1);
    Matrix<SIZE_N> vectorColumnCi = Matrix<SIZE_N>(N, 1);

    explicit Phase2(int n, int start = 1, int end = 100)
            : N(n)
    {
        genBi(start, end);
        genCi(start, end);
        genA1();
    }

    void genBi(int start, int end)
    {
        randGenColumnVector(vectorColumnBi, N, start, end);
    }

    void genCi(int start, int end)
    {
        randGenColumnVector(vectorColumnCi, N, start, end);
    }

    void genA1()
    {
        randGenMatrix(A1, N);
    }

    Matrix<SIZE_N> countY2()
    {
        Matrix<SIZE_N> result(N, 1);

        result = vectorColumnBi.sub(vectorColumnCi);
        result = A1.mul(result);

        return result;
    }
};


class Phase3
{
public:
    const int N;
    Matrix<SIZE_N> A2 = Matrix<SIZE_N>(N, N);
    Matrix<SIZE_N> B2 = Matrix<SIZE_N>(N, N);
    Matrix<SIZE_N> C2 = Matrix<SIZE_N>(N, N);

    explicit Phase3(int n, int start = 1, int end = 100)
            : N(n)
    {
        genA2(start, end);
        genB2(start, end);
        genC2();
    }

    void genA2(int start = 1, int end = 100)
    {
        randGenMatrix(A2, N, start, end);
    }

    void genB2(int start = 1, int end = 100)
    {
        randGenMatrix(B2, N, start, end);
    }

    void genC2()
    {
        for (int i = 0; i < N; ++i)
        {
            for (int j = 0; j < N; ++j)
            {
                C2[i][j] = 1 / (double)(i + j + 1);
            }
        }
    }

    Matrix<SIZE_N> countY3()
    {
        Matrix<SIZE_N> sub = B2.add(C2);
        return A2.mul(sub);
    }
};


int main()
{
    Phase1 phase1(4);
    Matrix<SIZE_N> Y1 = phase1.countY1();
    Matrix<SIZE_N> Y1_T = Y1.transpose();

    Phase2 phase2(4, 1, 100);
    Matrix<SIZE_N> Y2 = phase2.countY2();
    Matrix<SIZE_N> Y2_T = Y2.transpose();

    Phase3 phase3(4, 1, 100);
    Matrix<SIZE_N> Y3 = phase3.countY3();
    Matrix<SIZE_N> Y3_2 = Y3.mul(Y3);
    Matrix<SIZE_N> Y3_3 = Y3_2.mul(Y3);

    return 0;
}
