#include <vector>
#include <iostream>
#include <random>
#include <string>

#include <chrono>

typedef double type;

class Gauss
{
    std::vector<std::vector<type>> coefs;
    std::vector<std::vector<type>> coefs_tmp;
    std::vector<type> variables;
    std::vector<type> free_term;
    std::vector<type> free_term_tmp;
    int n, m;

public:
    Gauss(int _n, int _m, type **_coefs, type *_free_terms)
        : n(_n),
          m(_m),
          coefs(_m, std::vector<type>(_n)),
          coefs_tmp(_m, std::vector<type>(_n)),
          variables(_n),
          free_term(_m),
          free_term_tmp(_m)
    {
        for (int i = 0; i < _m; ++i)
        {
            for (int j = 0; j < _n; ++j)
            {
                coefs[i][j] = coefs_tmp[i][j] = _coefs[i][j];
            }
            free_term[i] = free_term_tmp[i] = _free_terms[i];
        }
    }
    double solve()
    {
        std::chrono::time_point start = std::chrono::high_resolution_clock::now();
        forward();
        backward();
        get_answer();
        std::chrono::time_point end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed_time = end - start;
        return elapsed_time.count();
    }

    void forward()
    {
        for (int i = 0; i < m; ++i)
        {
            for (int j = 0; j < n; ++j)
            {
                coefs_tmp[i][j] = coefs_tmp[i][j] / coefs[i][i];
            }

            free_term_tmp[i] = free_term_tmp[i] / coefs[i][i];

            for (int j = i + 1; j < n; ++j)
            {
                type coef = coefs_tmp[j][i] / coefs_tmp[i][i];
                for (int k = 0; k < n; ++k)
                {
                    coefs_tmp[j][k] = coefs_tmp[j][k] - coefs_tmp[i][k] * coef;
                }
                free_term_tmp[j] = free_term_tmp[j] - free_term_tmp[i] * coef;
            }

            for (int k = 0; k < m; ++k)
            {
                for (int j = 0; j < n; ++j)
                {
                    coefs[k][j] = coefs_tmp[k][j];
                }
                free_term[k] = free_term_tmp[k];
            }
        }
    }
    void backward()
    {
        for (int i = m - 1; i >= 0; --i)
        {
            for (int j = n - 1; j >= 0; --j)
            {
                coefs_tmp[i][j] = coefs_tmp[i][j] / coefs[i][i];
            }
            free_term_tmp[i] = free_term_tmp[i] / coefs[i][i];

            for (int j = i - 1; j >= 0; --j)
            {
                type coef = coefs_tmp[j][i] / coefs_tmp[i][i];
                for (int k = n - 1; k >= 0; --k)
                {
                    coefs_tmp[j][k] = coefs_tmp[j][k] - coefs_tmp[i][k] * coef;
                }
                free_term_tmp[j] = free_term_tmp[j] - free_term_tmp[i] * coef;
            }
        }
    }

    void get_answer()
    {
        for (int i = 0; i < n; ++i)
        {
            variables[i] = free_term_tmp[i];
        }
    }

    type &operator[](int ind)
    {
        return variables[ind];
    }
};

bool correctness_test_run()
{
    type line_1[] = {7, 1, 0};
    type line_2[] = {2, -10, 4};
    type line_3[] = {-10, -1, -8};

    type free_terms[] = {-47, -74, 148};

    type answer[] = {-7, 2, -10};

    type *coefs[] = {line_1, line_2, line_3};

    Gauss method{3, 3, coefs, free_terms};
    double test_time = method.solve();
    std::cout << "Elapsed sequential Gaussian time: " << test_time << " seconds\n";

    for (int i = 0; i < 3; ++i)
    {
        if (std::abs(answer[i] - method[i]) > 1e-9)
        {
            std::cout << "Small test failed. Calculations are incorrect.\n";
            return false;
        }
    }
    std::cout << "Small test passed.\n";
    return true;
}

bool performance_run(int task_size)
{
    type **coefs = new type *[task_size];
    type *free_term = new type[task_size];

    std::mt19937 generator{};

    for (int i = 0; i < task_size; ++i)
    {
        coefs[i] = new type[task_size];
        for (int j = 0; j < task_size; ++j)
        {
            coefs[i][j] = generator();
            if (i == j)
            {
                while (coefs[i][j] == 0)
                {
                    coefs[i][j] = generator();
                }
            }
        }
        free_term[i] = generator();
    }

    Gauss method(task_size, task_size, coefs, free_term);
    double seq_time = method.solve();
    std::cout << "Elapsed sequential Gaussian time: " << seq_time << " seconds\n";

    delete[] (free_term);
    for (int i = 0; i < task_size; ++i)
    {
        delete[] (coefs[i]);
    }
    delete[] (coefs);
    return true;
}

int main(int argc, char **argv)
{
    int task_size = 100;

    if (argc > 1)
    {
        task_size = std::stoi(argv[1]);
        std::cout << "Task size set to: " << task_size << std::endl;
    }
    else
    {
        std::cout << "Using default task size: " << task_size << std::endl;
    }
    correctness_test_run();
    performance_run(task_size);
}
