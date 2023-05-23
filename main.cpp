#include <random>
#include <iostream>
#include <iomanip>
#include <vector>
#include <map>

class SharedMatrix;

typedef long double CalVar;
typedef std::vector<CalVar> TVector;
typedef std::vector<TVector> TMatrix;
typedef std::vector<SharedMatrix> SMList;


class SharedMatrix
{
private:
    TMatrix _matrix;
    unsigned _rows;
    unsigned _cols;
    unsigned _flag;

public:
    SharedMatrix(unsigned rows, unsigned cols, unsigned flag)
    {
        this->_rows = rows;
        this->_cols = cols;
        this->_flag = flag;

        _matrix.resize(rows);
        for(int i = 0; i < rows; i++)
            _matrix[i].resize(cols);
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

    unsigned flag()
    {
        return _flag;
    }
};


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

    bool compare(Matrix &oth)
    {
        if (_rows != oth.rows() || _cols != oth.cols())
            return false;

        auto &othMatrix = oth.matrix();

        for(int i = 0; i < _rows; i++)
            for(int j = 0; j < _cols; j++)
                if(_matrix[i][j] != othMatrix[i][j])
                    return false;

        return true;
    }

    SMList splitIntoMatricesRow(int parts)
    {
        std::vector<SharedMatrix> result;

        int rowsPerPart = _rows / parts;
        int rowsLeft = _rows % parts;

        int start = 0;
        int end = rowsPerPart;

        for(int i = 0; i < parts; i++)
        {
            if (rowsLeft > 0)
            {
                end++;
                rowsLeft--;
            }

            result.push_back(SharedMatrix(end - start, _cols, i));
            start = end;
            end += rowsPerPart;
        }

        return result;
    }

    SMList splitIntoMatricesCol(int parts)
    {
        std::vector<SharedMatrix> result;

        int colsPerPart = _cols / parts;
        int colsLeft = _cols % parts;

        int start = 0;
        int end = colsPerPart;

        for(int i = 0; i < parts; i++)
        {
            if (colsLeft > 0)
            {
                end++;
                colsLeft--;
            }

            result.push_back(SharedMatrix(_rows, end - start, i));
            start = end;
            end += colsPerPart;
        }

        return result;
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

    for(int i = 0; i < matrixA.rows()-50; i++)
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
    Matrix matrixA(60, 248);
    Matrix matrixB(248, 149);

    genMatrix(matrixA);
    genMatrix(matrixB);

    // 7 -> 4 -> 8 -> 6 -> 2 -> 1 -> 3 -> 5
    std::map<int, int> communicationData {
        {7, 4},
        {4, 8},
        {8, 6},
        {6, 2},
        {2, 1},
        {1, 3},
        {3, 5},
        {5, 7}
    };

    auto vec = matrixA.splitIntoMatricesRow(8);

    for (auto &el: vec)
    {
        std::cout << el.flag() << " " << el.rows() << "|" << el.cols() << std::endl;
    }

    std::cout << std::endl;

    auto vec2 = matrixB.splitIntoMatricesCol(8);
    for (auto &el: vec2)
    {
        std::cout << el.flag() << " " << el.rows() << "|" << el.cols() << std::endl;
    }

    return 0;
}