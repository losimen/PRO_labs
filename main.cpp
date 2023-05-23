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
    {
        for(int j = 0; j < matrix.cols(); j++)
        {
            m_matrix[i][j] = randGenNumber();
        }
    }

}


int main()
{
    Matrix matrixA(6, 3);
    Matrix matrixB(5, 6);

    genMatrix(matrixA);
    matrixA.print();

    return 0;
}