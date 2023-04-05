#include <vector>
#include <random>
#include <iostream>

class Phase1
{
public:
    const int N;
    std::vector<int> vector_column_bi {};
    std::vector<std::vector<int>> A {};

    explicit Phase1(int n)
            : N(n)
    {}

    std::vector<int> countY1()
    {
        std::vector<int> result;
        result.reserve(vector_column_bi.size());

        for (const auto &row: A)
        {
            int sum = 0;
            for (unsigned i = 0; i < row.size(); ++i)
            {
                sum += row[i] * vector_column_bi[i];
            }
            result.emplace_back(sum);
        }

        return result;
    }

    void genBi(int num)
    {
        vector_column_bi.reserve(N);

        for (unsigned i = 0; i < N; ++i)
        {
            vector_column_bi.emplace_back(int(num * (i + 1)));
        }
    }

    void genA()
    {
        std::random_device rd;
        std::mt19937 mt(rd());
        std::uniform_real_distribution<double> dist(1, 100);

        A.reserve(N);
        for (unsigned i = 0; i < N; ++i)
        {
            std::vector<int> row;
            row.reserve(N);

            for (unsigned j = 0; j < N; ++j)
            {
                row.emplace_back((int)dist(mt));
            }

            A.emplace_back(row);
        }
    }
};


class Phase2
{
public:
    const int N;
    std::vector<std::vector<int>> A1;
    std::vector<int> vector_column_bi;
    std::vector<int> vector_column_ci;

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
        std::random_device rd;
        std::mt19937 mt(rd());
        std::uniform_real_distribution<double> dist(1, 100);

        A1.reserve(N);
        for (unsigned i = 0; i < N; ++i)
        {
            std::vector<int> row;
            row.reserve(N);

            for (unsigned j = 0; j < N; ++j)
            {
                row.emplace_back((int)dist(mt));
            }

            A1.emplace_back(row);
        }
    }

    std::vector<int> countY2()
    {
        auto subTwoVectors = [](const std::vector<int> &v1, const std::vector<int> &v2)
        {
            std::vector<int> result;
            result.reserve(v1.size());

            for (unsigned i = 0; i < v1.size(); ++i)
            {
                result.emplace_back(v1[i] - v2[i]);
            }

            return result;
        };

        std::vector<int> result;
        result.reserve(vector_column_ci.size());

        for (const auto &row: A1)
        {
            const auto sub = subTwoVectors(row, vector_column_bi);

            int sum = 0;
            for (unsigned i = 0; i < sub.size(); ++i)
            {
                sum += sub[i] * vector_column_ci[i];
            }

            result.emplace_back(sum);
        }

        return result;
    }
};







int main()
{
    const int N = 4;

    Phase1 p1(N);
    p1.genA();
    p1.genBi(9);
    std::vector<int> y1 = p1.countY1();

    for (const auto &elem: y1)
        std::cout << elem << " ";
    std::cout << std::endl;

    Phase2 p2(N);
    p2.genA1();
    p2.genBi(9);
    p2.genCi(7);

    std::vector<int> y2 = p2.countY2();
    for (const auto &elem: y2)
        std::cout << elem << " ";
    std::cout << std::endl;

//    for (const auto &row: p1.A)
//    {
//        for (const auto &elem: row)
//        {
//            std::cout << elem << " ";
//        }
//        std::cout << std::endl;
//    }

    return 0;
}
