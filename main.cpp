#include <random>
#include <iostream>
#include <iomanip>
#include <vector>

typedef long double CalVar;
typedef std::vector<CalVar> TVector;
typedef std::vector<TVector> TMatrix;


class Matrix
{
private:
    TMatrix _matrix;
    unsigned _rows;
    unsigned _cols;

public:
    Matrix(unsigned rows, unsigned cols)
    {
        this->_rows = rows;
        this->_cols = cols;

        _matrix.resize(rows);
        for(int i = 0; i < rows; i++)
        {
            _matrix[i].resize(cols);
        }
    }

    TMatrix &matrix()
    {
        return _matrix;
    }

    unsigned rows()
    {
        return _rows;
    }

    unsigned cols()
    {
        return _cols;
    }

    void print()
    {
        for(int i = 0; i < _rows; i++)
        {
            for(int j = 0; j < _cols; j++)
            {
                std::cout << std::setw(5) << _matrix[i][j] << " ";
            }
            std::cout << std::endl;
        }
    }
};


CalVar randGenNumber(int start = 1, int end = 100)
{
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<CalVar> dist(start, end);

    return dist(mt);
}


void genMatrix(Matrix &matrix)
{
    auto &m_matrix = matrix.matrix();

    for(int i = 0; i < matrix.rows(); i++)
        for(int j = 0; j < matrix.cols(); j++)
            m_matrix[i][j] = randGenNumber();
}


void mulMatrices(Matrix &result, Matrix &matrixA, Matrix &matrixB)
{
    auto m_matrixA = matrixA.matrix();
    auto m_matrixB = matrixB.matrix();
    auto &m_result = result.matrix();

    for(int i = 0; i < matrixA.rows(); i++)
    {
        for(int j = 0; j < matrixB.cols(); j++)
        {
            m_result[i][j] = 0;
            for(int k = 0; k < matrixA.cols(); k++)
            {
                m_result[i][j] += m_matrixA[i][k] * m_matrixB[k][j];
            }
        }
    }
}


int main()
{
    Matrix matrixA(10, 5);
    Matrix matrixB(5, 10);

    genMatrix(matrixA);
    genMatrix(matrixB);

//
//    std::cout << "Matrix A:" << std::endl;
//    matrixA.print();
//    std::cout << std::endl;
//
//    std::cout << "Matrix B:" << std::endl;
//    matrixB.print();
//    std::cout << std::endl;

    Matrix result(matrixA.rows(), matrixB.cols());
    mulMatrices(result, matrixA, matrixB);
    std::cout << "Result:" << std::endl;
    result.print();
    std::cout << std::endl;

    return 0;
}