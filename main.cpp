#include <random>
#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <mpi.h>
#include <thread>
#include <chrono>
#include <fstream>

class SharedMatrix;

typedef long double CalVar;
typedef std::vector<CalVar> TVector;
typedef std::vector<TVector> TMatrix;
typedef std::vector<SharedMatrix> SMVec;


// 7 -> 4 -> 8 -> 6 -> 2 -> 1 -> 3 -> 5
std::map<int, int> communicationMulSend {
        {7, 4},
        {4, 8},
        {8, 6},
        {6, 2},
        {2, 1},
        {1, 3},
        {3, 5},
        {5, 7}
};


std::map<int, int> communicationMulRecv {
        {7, 5},
        {4, 7},
        {8, 4},
        {6, 8},
        {2, 6},
        {1, 2},
        {3, 1},
        {5, 3}
};


std::map<int, int> communicationDataSend {
        {7, 4},
        {4, 8},
        {8, 6},
        {6, 2},
        {2, 1},
        {1, 3},
        {3, 5},
        {5, 0},
        {0, 7}
};

std::map<int, int> communicationDataRecv {
        {0, 5},
        {7, 0},
        {4, 7},
        {8, 4},
        {6, 8},
        {2, 6},
        {1, 2},
        {3, 1},
        {5, 3}
};





class SharedMatrix
{
private:
    TMatrix _matrix;
    unsigned _rows;
    unsigned _cols;
    unsigned _flag;
    unsigned _statusCode = 0;

public:
    SharedMatrix()
    {
        this->_rows = 0;
        this->_cols = 0;
        this->_flag = 0;
    }

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

    unsigned statusCode()
    {
        return _statusCode;
    }

    void setStatusCode(unsigned statusCode)
    {
        this->_statusCode = statusCode;
    }

    CalVar *serializeData()
    {
        // status + flag + rows + cols + matrix
        CalVar *result = new CalVar[4+(_rows * _cols)];

        result[0] = _statusCode;
        result[1] = _flag;
        result[2] = _rows;
        result[3] = _cols;

        int index = 4;
        for(int i = 0; i < _rows; i++)
            for(int j = 0; j < _cols; j++)
                result[index++] = _matrix[i][j];

        return result;
    }

    void parseData(CalVar *data)
    {
        _statusCode = data[0];
        _flag = data[1];
        _rows = data[2];
        _cols = data[3];

        _matrix.resize(_rows);
        for(int i = 0; i < _rows; i++)
            _matrix[i].resize(_cols);

        int index = 4;
        for(int i = 0; i < _rows; i++)
            for(int j = 0; j < _cols; j++)
                _matrix[i][j] = data[index++];
    }

    void print()
    {
        for (int i = 0; i < _rows; i++)
        {
            for (int j = 0; j < _cols; j++)
                std::cout << std::setw(5) << _matrix[i][j] << " ";
            std::cout << std::endl;
        }
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

    SMVec splitIntoMatricesRow(int parts)
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

            auto newM = SharedMatrix(end - start, _cols, i);
            newM.matrix() = TMatrix(_matrix.begin() + start, _matrix.begin() + end);
            result.push_back(newM);

            start = end;
            end += rowsPerPart;
        }

        return result;
    }

    SMVec splitIntoMatricesCol(int parts)
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

            auto newM = SharedMatrix(_rows, end - start, i);
            for(int j = 0; j < _rows; j++)
                newM.matrix()[j] = TVector(_matrix[j].begin() + start, _matrix[j].begin() + end);
            result.push_back(newM);

            start = end;
            end += colsPerPart;
        }

        return result;
    }

    void insertBlock(SharedMatrix &block, int row, int col)
    {
        auto &blockMatrix = block.matrix();

        for(int i = 0; i < block.rows(); i++)
            for(int j = 0; j < block.cols(); j++)
                _matrix[row + i][col + j] = blockMatrix[i][j];
    }

    void writeToFile(std::string fileName)
    {
        std::ofstream file(fileName);

        file << _rows << " " << _cols << std::endl;

        for(int i = 0; i < _rows; i++)
        {
            for(int j = 0; j < _cols; j++)
                file << _matrix[i][j] << " ";
            file << std::endl;
        }

        file.close();
    }
};


class ProcessCommunicator
{
public:
    static void send(int rankToSend, SharedMatrix &matrix)
    {
        auto data = matrix.serializeData();
        const int size = 4 + (matrix.rows() * matrix.cols());

        MPI_Send(data, size, MPI_LONG_DOUBLE, rankToSend, 0, MPI_COMM_WORLD);

        delete[] data;
    }

    static SharedMatrix recv(int recvFrom = MPI_ANY_SOURCE)
    {
        MPI_Status status;
        MPI_Probe(recvFrom, 0, MPI_COMM_WORLD, &status);

        int size;
        MPI_Get_count(&status, MPI_LONG_DOUBLE, &size);

        CalVar *data = new CalVar[size];
        MPI_Recv(data, size, MPI_LONG_DOUBLE, recvFrom, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        SharedMatrix matrix;
        matrix.parseData(data);

        delete[] data;

        return matrix;
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


void mulSharedMatrices(SharedMatrix &result, SharedMatrix &matrixA, SharedMatrix &matrixB)
{
    auto &m_matrixA = matrixA.matrix();
    auto &m_matrixB = matrixB.matrix();
    auto &m_result = result.matrix();

    if (matrixA.cols() != matrixB.rows())
        throw std::runtime_error("Matrix A cols != Matrix B rows");

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


void jobRank0(int procRank)
{
    Matrix matrixA(60, 248);
    Matrix matrixB(248, 149);

    Matrix resultP(60, 149);
    Matrix resultS(60, 149);

    SharedMatrix zero(0, 0, 15);
    zero.setStatusCode(15);

    genMatrix(matrixA);
    genMatrix(matrixB);

    auto vecSM_A = matrixA.splitIntoMatricesRow(8);
    auto vecSM_B = matrixB.splitIntoMatricesCol(8);

    for(unsigned i = 0; i < vecSM_A.size(); ++i)
    {
        auto el = vecSM_A[i];

        el.setStatusCode(i == vecSM_A.size() - 1 ? 1 : 0);

        ProcessCommunicator::send(communicationDataSend[procRank], el);
        ProcessCommunicator::recv(communicationDataRecv[procRank]);
    }

    for (unsigned i = 0; i < vecSM_B.size(); ++i)
    {
        auto el = vecSM_B[i];

        el.setStatusCode(i == vecSM_B.size() - 1 ? 1 : 0);

        ProcessCommunicator::send(communicationDataSend[procRank], el);
        ProcessCommunicator::recv(communicationDataRecv[procRank]);
    }

    for (unsigned i = 0; i < 8; ++i)
    {
        for (unsigned j = 1; j <= 8; ++j)
        {
            auto matrix = ProcessCommunicator::recv(j);

            unsigned a = matrix.flag();
            unsigned b = matrix.statusCode();

            resultP.insertBlock(matrix, a*7, b*7);

            for (unsigned k = 1; k <= 8; ++k)
            {
                ProcessCommunicator::send(k, zero);
            }
        }
    }

    std::cout << "Result: " << std::endl;
    mulMatrices(resultS, matrixA, matrixB);

    resultP.writeToFile("resultP.txt");
    resultS.writeToFile("resultS.txt");
}


void jobRankN(int procRank)
{
    unsigned status = 0;
    SharedMatrix ma;
    SharedMatrix mb;

    // recv ma;
    while (status == 0)
    {
        auto matrix = ProcessCommunicator::recv(communicationDataRecv[procRank]);

        if (matrix.flag() == procRank-1)
        {
            ma = matrix;
        }

        ProcessCommunicator::send(communicationDataSend[procRank], matrix);
        status = matrix.statusCode();
    }

    status = 0;
    while (status == 0)
    {
        auto matrix = ProcessCommunicator::recv(communicationDataRecv[procRank]);

        if (matrix.flag() == procRank-1)
        {
            mb = matrix;
        }

        ProcessCommunicator::send(communicationDataSend[procRank], matrix);
        status = matrix.statusCode();
    }

    for (unsigned i = 0; i < 8; ++i)
    {
        auto result = SharedMatrix(ma.rows(), mb.cols(), procRank-1);
        result.setStatusCode(mb.flag());
        mulSharedMatrices(result, ma, mb);

//        std::cout << procRank << "-result " << result.flag() << " " << result.rows() << "|" << result.cols() << std::endl;
        ProcessCommunicator::send(0, result);
        auto zero = ProcessCommunicator::recv(0);
//        std::cout << procRank << "-sent " << zero.flag() << " " << zero.rows() << "|" << zero.cols() << std::endl;

        if (procRank == 8)
        {
            mb = ProcessCommunicator::recv(communicationMulRecv[procRank]);
            ProcessCommunicator::send(communicationMulSend[procRank], mb);
        }
        else
        {
            ProcessCommunicator::send(communicationMulSend[procRank], mb);
            mb = ProcessCommunicator::recv(communicationMulRecv[procRank]);
        }
    }
}


int main(int argc, char* argv[])
{
    int procRank;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);

    if (procRank == 0)
    {
        jobRank0(procRank);
    }
    else
    {
        jobRankN(procRank);
    }

    MPI_Finalize();
    return 0;
}