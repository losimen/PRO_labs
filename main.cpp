#include <random>
#include <iostream>
#include <array>
#include <iomanip>
#include <map>
#include <mpi.h>
#include <fstream>

constexpr unsigned SIZE_N = 4;

typedef long double CalVar;
typedef std::array<std::array<CalVar, SIZE_N>, SIZE_N> TMatrix;
typedef std::array<CalVar, SIZE_N> TVector;


template <unsigned N>
class Matrix
{
private:
    TMatrix m_matrix{};

    Matrix mul(Matrix &oth)
    {
        if (this->cols != oth.rows)
            throw std::invalid_argument("m_matrix sizes are not equal");

        Matrix<SIZE_N> result(this->rows, oth.cols);

        for (size_t i = 0; i < rows; i++) {
            for (size_t j = 0; j < oth.cols; j++) {
                for (size_t k = 0; k < cols; k++) {
                    result[i][j] += m_matrix[i][k] * oth.m_matrix[k][j];
                }
            }
        }

        return result;
    }

    Matrix mul(CalVar oth)
    {
        Matrix<SIZE_N> result(this->rows, this->cols);

        for (size_t i = 0; i < rows; i++) {
            for (size_t j = 0; j < cols; j++) {
                result[i][j] = m_matrix[i][j] * oth;
            }
        }

        return result;
    }

    Matrix sub(Matrix &oth)
    {
        if (this->rows != oth.rows && this->cols != oth.cols)
            throw std::invalid_argument("m_matrix sizes are not equal");

        Matrix<SIZE_N> result(std::max(this->rows, oth.rows), std::max(this->cols, oth.cols));

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
        if (this->rows != oth.rows && this->cols != oth.cols)
            throw std::invalid_argument("m_matrix sizes are not equal");

        Matrix<SIZE_N> result(std::max(this->rows, oth.rows), std::max(this->cols, oth.cols));

        for (unsigned i = 0; i < m_matrix.size(); ++i)
        {
            for (unsigned j = 0; j < m_matrix[i].size(); ++j)
            {
                result[i][j] = m_matrix[i][j] + oth[i][j];
            }
        }

        return result;
    }
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

    Matrix transpose()
    {
        Matrix<SIZE_N> result(this->cols, this->rows);

        for (unsigned i = 0; i < m_matrix.size(); ++i)
            for (unsigned j = 0; j < m_matrix[i].size(); ++j)
                result[i][j] = m_matrix[j][i];

        return result;
    }

    template<typename T>
    Matrix operator*(T &&oth)
    {
        return mul(oth);
    }

    template<typename T>
    Matrix operator-(T &&oth)
    {
        return sub(oth);
    }

    template<typename T>
    Matrix operator+(T &&oth)
    {
        return add(oth);
    }

    void print()
    {
        std::cout << std::fixed;
        for (auto & i : m_matrix)
        {
            for (long double j : i)
            {
                std::cout << std::setw(12) << std::setprecision(3) << j << " ";
            }
            std::cout << std::endl;
        }
    }

    void writeToFile(const std::string &fileName)
    {
        std::ofstream file(fileName);
        file << std::fixed;
        for (auto & i : m_matrix)
        {
            for (long double j : i)
            {
                file << j << " ";
            }
            file << std::endl;
        }
    }

    void readFromFile(const std::string &fileName)
    {
        std::ifstream file(fileName);
        for (auto & i : m_matrix)
        {
            for (long double & j : i)
            {
                file >> j;
            }
        }
    }
};

void randGenMatrix(Matrix<SIZE_N> &A, int N, int start = 1, int end = 100)
{
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<CalVar> dist(start, end);

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
    std::uniform_real_distribution<CalVar> dist(start, end);

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

    explicit Phase1(int n, int start = 1, int end = 100)
            : N(n)
    {
        genBi(9);
        genA(start, end);
    }

    Matrix<SIZE_N> countY1()
    {
        Matrix<SIZE_N> result(N, 1);

        for (unsigned i = 0; i < N; ++i)
        {
            CalVar sum = 0;
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

        result = vectorColumnBi - vectorColumnCi;
        result = A1 * result;

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
                C2[i][j] = 1 / (CalVar)(i + j + 1);
            }
        }
    }

    Matrix<SIZE_N> countY3()
    {
        Matrix<SIZE_N> sub = B2 + C2;
        return A2 * sub;
    }
};


class Process
{
private:
    static void serializeMatrix(Matrix<SIZE_N>& matrix, std::vector<CalVar>& vec)
    {
        for (int i = 0; i < SIZE_N; ++i)
        {
            for (int j = 0; j < SIZE_N; ++j)
            {
                vec[i*SIZE_N + j] = matrix[i][j];
            }
        }
    }

    static void deserializeMatrix(Matrix<SIZE_N>& matrix, std::array<CalVar, SIZE_N*SIZE_N>& vec)
    {
        for (int i = 0; i < SIZE_N; ++i)
        {
            for (int j = 0; j < SIZE_N; ++j)
            {
                matrix[i][j] = vec[i * SIZE_N + j];
            }
        }
    }


public:
    static void sendMatrix(int rankToSend, Matrix<SIZE_N>& matrix)
    {
        std::vector<CalVar> vecToSend(SIZE_N * SIZE_N);
        Process::serializeMatrix(matrix, vecToSend);

        MPI_Send(vecToSend.data(), SIZE_N * SIZE_N, MPI_LONG_DOUBLE, rankToSend, 0,
                 MPI_COMM_WORLD);
    }

    static Matrix<SIZE_N> recvMatrix(int recvFrom)
    {
        Matrix<SIZE_N> matrix(4, 4);
        std::array<CalVar , SIZE_N*SIZE_N> recvMatrix{};
        MPI_Status status;

        MPI_Recv(recvMatrix.data(), SIZE_N * SIZE_N, MPI_LONG_DOUBLE,
                 recvFrom, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

        deserializeMatrix(matrix, recvMatrix);
        return matrix;
    }

    static void sendVector(int rankToSend, std::array<CalVar, SIZE_N>& vec)
    {
        MPI_Send(vec.data(), SIZE_N, MPI_LONG_DOUBLE, rankToSend, 0,
                 MPI_COMM_WORLD);
    }

    static Matrix<SIZE_N> recvVector(int recvFrom)
    {
        std::array<CalVar , SIZE_N> rec{};
        MPI_Status status;

        MPI_Recv(rec.data(), SIZE_N, MPI_LONG_DOUBLE,
                 recvFrom, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

        Matrix<SIZE_N> matrix(1, 4);
        for (int i = 0; i < SIZE_N; ++i)
        {
            matrix[0][i] = rec[i];
        }

        return matrix;
    }

    static void sendScalar(int rankToSend, CalVar scalar)
    {
        MPI_Send(&scalar, 1, MPI_LONG_DOUBLE, rankToSend, 0,
                 MPI_COMM_WORLD);
    }

    static CalVar recvScalar(int recvFrom)
    {
        CalVar scalar;
        MPI_Status status;

        MPI_Recv(&scalar, 1, MPI_LONG_DOUBLE,
                 recvFrom, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

        return scalar;
    }
};


Matrix<SIZE_N> FirstBraces(Matrix<SIZE_N> &Y3, Matrix<SIZE_N> &Y3_2,
                           Matrix<SIZE_N> &Y1, Matrix<SIZE_N> &Y1_T,
                           Matrix<SIZE_N> &Y2, Matrix<SIZE_N> &Y2_T,
                           CalVar a)
{
    Matrix<SIZE_N> result = Y3 + Y3_2 * a + Y2 * Y2_T;
    return result;
}

Matrix<SIZE_N> SecondBraces(Matrix<SIZE_N> &Y3, Matrix<SIZE_N> &Y3_3, CalVar a)
{
    Matrix<SIZE_N> result = Y3 * a + Y3_3;
    return result;
}


void countEquation(Phase1 &phase1, Phase2 &phase2, Phase3 &phase3)
{
    Matrix<SIZE_N> Y1 = phase1.countY1();
    Matrix<SIZE_N> Y1_T = Y1.transpose();
    CalVar a = (Y1_T * Y1)[0][0];

    Matrix<SIZE_N> Y2 = phase2.countY2();
    Matrix<SIZE_N> Y2_T = Y2.transpose();

    Matrix<SIZE_N> Y3 = phase3.countY3();
    Matrix<SIZE_N> Y3_2 = Y3 * Y3;
    Matrix<SIZE_N> Y3_3 = Y3_2 * Y3;

    Matrix<SIZE_N> firstBraces = FirstBraces(Y3, Y3_2, Y1, Y1_T, Y2, Y2_T, a);
    Matrix<SIZE_N> secondBraces = SecondBraces(Y3, Y3_3, a);

    secondBraces = Y1_T * secondBraces;
    firstBraces = Y2_T * firstBraces;

    auto final = secondBraces + firstBraces;
    std::cout << "Final result one process: " << std::endl;
    final.print();
}


void doFirstProcess(Phase1 &phase1)
{
    Matrix<SIZE_N> Y1 = phase1.countY1();
    Matrix<SIZE_N> Y1_T = Y1.transpose();
    CalVar a = (Y1_T * Y1)[0][0];

    Process::sendScalar(1, a);
    Matrix<SIZE_N> Y3 = Process::recvMatrix(2);
    Matrix<SIZE_N> Y3_3 = Process::recvMatrix(2);
    Matrix<SIZE_N> Y3_Mul_a = Y3 * a + Y3_3;

    Matrix<SIZE_N> Y1_T_Mul_Y3_Mul_a_Plus_Y3_3 = Y1_T * Y3_Mul_a;

    Process::sendVector(1, Y1_T_Mul_Y3_Mul_a_Plus_Y3_3[0]);
}

void doSecondProcess(Phase2 &phase2)
{
    Matrix<SIZE_N> Y2 = phase2.countY2();
    Matrix<SIZE_N> Y2_T = Y2.transpose();

    // matrix
    Matrix<SIZE_N> Y2_Mul_Y2_T = Y2 * Y2_T;

    CalVar a = Process::recvScalar(0);
    Matrix<SIZE_N> Y3_2 = Process::recvMatrix(2);

    Matrix<SIZE_N> Y3_2_Mul_a = Y3_2 * a;
    Matrix<SIZE_N> Y3_2_Mul_a_Plus_Y2_Mul_Y2_T = Y3_2_Mul_a + Y2_Mul_Y2_T;

    Matrix<SIZE_N> Y3 = Process::recvMatrix(2);
    Matrix<SIZE_N> Y3_2_Mul_a_Plus_Y2_Mul_Y2_T_Plus_Y3 = Y3_2_Mul_a_Plus_Y2_Mul_Y2_T + Y3;
    Matrix<SIZE_N> fin = Y2_T * Y3_2_Mul_a_Plus_Y2_Mul_Y2_T_Plus_Y3;

    Matrix<SIZE_N> recived = Process::recvVector(0);

    std::cout << "Final result two process: " << std::endl;
    (fin + recived).print();
}

void doThirdProcess(Phase3 &phase3)
{
    Matrix<SIZE_N> Y3 = phase3.countY3();
    Process::sendMatrix(0, Y3);

    Matrix<SIZE_N> Y3_2 = Y3 * Y3;
    Process::sendMatrix(1, Y3_2);

    Matrix<SIZE_N> Y3_3 = Y3_2 * Y3;
    Process::sendMatrix(0, Y3_3);
    Process::sendMatrix(1, Y3);
}


void saveToFile()
{
    Phase1 phase1(4, 1, 3);
    Phase2 phase2(4, 1, 3);
    Phase3 phase3(4, 1, 3);

    phase1.vectorColumnBi.writeToFile("vectorColumnBi_1.txt");
    phase1.A.writeToFile("A.txt");

    phase2.A1.writeToFile("A1.txt");
    phase2.vectorColumnBi.writeToFile("vectorColumnBi_2.txt");
    phase2.vectorColumnCi.writeToFile("vectorColumnCi.txt");

    phase3.A2.writeToFile("A2.txt");
    phase3.B2.writeToFile("B2.txt");
    phase3.C2.writeToFile("C2.txt");
}


void readFromFile(Phase1 &phase1, Phase2 &phase2, Phase3 &phase3)
{
    phase1.vectorColumnBi.readFromFile("vectorColumnBi_1.txt");
    phase1.vectorColumnBi.rows = 4;
    phase1.vectorColumnBi.cols = 1;

    phase1.A.readFromFile("A.txt");
    phase1.A.rows = 4;
    phase1.A.cols = 4;

    phase2.A1.readFromFile("A1.txt");
    phase2.A1.rows = 4;
    phase2.A1.cols = 4;

    phase2.vectorColumnBi.readFromFile("vectorColumnBi_2.txt");
    phase2.vectorColumnBi.rows = 4;
    phase2.vectorColumnBi.cols = 1;

    phase2.vectorColumnCi.readFromFile("vectorColumnCi.txt");
    phase2.vectorColumnCi.rows = 4;
    phase2.vectorColumnCi.cols = 1;

    phase3.A2.readFromFile("A2.txt");
    phase3.A2.rows = 4;
    phase3.A2.cols = 4;

    phase3.B2.readFromFile("B2.txt");
    phase3.B2.rows = 4;
    phase3.B2.cols = 4;

    phase3.C2.readFromFile("C2.txt");
    phase3.C2.rows = 4;
    phase3.C2.cols = 4;
}


int main(int argc, char* argv[])
{
    int procRank;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    Phase1 phase1(4, 1, 3);
    Phase2 phase2(4, 1, 3);
    Phase3 phase3(4, 1, 3);

    readFromFile(phase1, phase2, phase3);

    if (procRank == 0)
        doFirstProcess(phase1);
    else if (procRank == 1)
        doSecondProcess(phase2);
    else if (procRank == 2)
        doThirdProcess(phase3);
    else
        countEquation(phase1, phase2, phase3);

    MPI_Finalize();
    return 0;
}
