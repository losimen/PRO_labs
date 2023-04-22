#include <random>
#include <iostream>
#include <array>
#include <iomanip>
#include <map>
#include <fstream>

constexpr unsigned SIZE_N = 6;

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
                std::cout << std::setw(4) << std::setprecision(0) << j << " ";
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

void genMatrixA(Matrix<SIZE_N> &A)
{
    for(int i = 0; i < SIZE_N; i++) {
        for(int j = 0; j < SIZE_N; j++) {
            if(i < SIZE_N / 2)
            {
                if(j < i || j >= SIZE_N - i)
                    A[i][j] = 0;
                else
                    A[i][j] = 1;
            }
            else {
                if(j < SIZE_N - i - 1 || j > i)
                    A[i][j] = 0;
                else
                    A[i][j] = 1;
            }
        }
    }
}

int main()
{
    Matrix<SIZE_N> A(SIZE_N, SIZE_N);
    genMatrixA(A);
    A.print();

    return 0;
}