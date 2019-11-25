#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <cmath>
#include <future>
#include <iostream>
#include <random>
#include <string>
#include <vector>


using matrix = std::vector<std::vector<double>>;
const double EPSILON = 1e-16;
const double EPSILON_RESIDUAL = 1e-10;
const int DIMENSION = 1000;
const std::pair<int, int> RANDOM_RANGE = {1, 10000};

class ZeroDiagonalEception : std::exception {};

void calculate_column(int k, matrix &A, matrix &X);
void decompose_to_LU(matrix &a_matrix);
matrix create_identity_matrix(int dim);
matrix inverse_matrix(matrix &mtx);
void print_matrix(const matrix &mtx);
matrix create_random_matrix(int dim);
matrix transpose_matrix(const matrix &mtx);
matrix substract_matrices(const matrix &a_matrix, const matrix &b_matrix);
matrix multiply_matrices(const matrix &a_matrix, const matrix &b_matrix);
matrix multiply_by_scalar(const matrix &a_matrix, double val);
double calculate_matrix_trace(const matrix &mtx);
matrix calculate_R_matrix(const matrix &I_matrix, const matrix &BxA_matrix);
bool check_R_matrix(const matrix &R_matrix);
matrix calculate_first_B_matrix(const matrix &trans_A_matrix,
                                const matrix &A_matrix);
matrix calculate_next_B_matrix(const matrix &B_matrix,
                               const matrix &BxA_matrix);
matrix inverse_matrix_iterative(const matrix &A_matrix);

/*
 * Prints out the matrix.
 * Use to debug.
 */
void print_matrix(const matrix &mtx) {
    for (const auto &raw : mtx) {
        for (const auto &element : raw)
            std::cout << element << " ";
        std::cout << "\n";
    }
}

/*
 * Create random square matrix. Its dimension is passed in argument.
 * Random values range is hardcoded in RANDOM_RANGE const variable.
 */
matrix create_random_matrix(int dim) {
    std::default_random_engine generator;
    std::uniform_int_distribution<int> dis(RANDOM_RANGE.first,
                                           RANDOM_RANGE.second);
    matrix random_matrix((dim), std::vector<double>(dim));
    for (auto &raw : random_matrix)
        for (auto &element : raw)
            element = dis(generator);

    return random_matrix;
}

/*
 * Create matrix with ones on the main diagonal and zeros elsewhere.
 */
matrix create_identity_matrix(int dim) {
    matrix i_matrix((dim), std::vector<double>(dim));
    i_matrix[0][0] = 1.0;
    for (int n = 0; n < dim; ++n)
        for (int m = 0; m < dim; ++m)
            if (n == m)
                i_matrix[n][m] = 1.0;
    return i_matrix;
}

/*
 * Factorization of a given square matrix into two triangular matrices L and U.
 * Matrices are merged in the given square matrix.
 */
void decompose_to_LU(matrix &a_matrix) {
    for (int k = 0; k < a_matrix.size() - 1; ++k) {
        if (fabs(a_matrix[k][k]) < EPSILON)
            throw ZeroDiagonalEception();

        for (int i = k + 1; i < a_matrix[0].size(); ++i)
            a_matrix[i][k] /= a_matrix[k][k];

        for (int i = k + 1; i < a_matrix.size(); ++i)
            for (int j = k + 1; j < a_matrix[0].size(); ++j)
                a_matrix[i][j] -= a_matrix[i][k] * a_matrix[k][j];
    }
}

/*
 * Calculate column of inversed matrix.
 */
void calculate_column(int k, matrix &a_matrix, matrix &x_matrix) {
    for (int i = 1; i < a_matrix.size(); ++i) {
        double s = 0;

        for (int j = 0; j < i; ++j)
            s += a_matrix[i][j] * x_matrix[j][k];

        x_matrix[i][k] -= s;
    }

    if (fabs(a_matrix[a_matrix.size() - 1][a_matrix[0].size() - 1]) < EPSILON)
        throw ZeroDiagonalEception();

    x_matrix[x_matrix.size() - 1][k] /=
        a_matrix[a_matrix.size() - 1][a_matrix[0].size() - 1];

    for (int i = x_matrix.size() - 2; i >= 0; --i) {
        double s = 0;

        for (int j = i + 1; j < a_matrix[0].size(); ++j)
            s += a_matrix[i][j] * x_matrix[j][k];

        if (fabs(a_matrix[i][i]) < EPSILON)
            throw ZeroDiagonalEception();

        x_matrix[i][k] = (x_matrix[i][k] - s) / a_matrix[i][i];
    }
}

/*
 * Inverses given matrix.
 */
matrix inverse_matrix(matrix &a_matrix) {
    matrix i_matrix = create_identity_matrix(a_matrix.size());
    try {
        decompose_to_LU(a_matrix);
        for (int i = 0; i < a_matrix.size(); ++i)
            calculate_column(i, a_matrix, i_matrix);
    } catch (ZeroDiagonalEception &) {
        std::cout << "Division by zero exception!\n"
                  << "There was a zero on the diagonal of a matrix.\n";
        exit(EXIT_FAILURE);
    }

    return i_matrix;
}

/*
 * Transposes square matrix.
 */
matrix transpose_matrix(const matrix &mtx) {
    matrix t_matrix((mtx.size()), std::vector<double>(mtx.size()));
    for (size_t i = 0; i < mtx.size(); ++i)
        for (size_t j = 0; j < mtx[0].size(); ++j)
            t_matrix[j][i] = mtx[i][j];
    return t_matrix;
}

/*
 * Substracts matrix A from matrix B
 */
matrix substract_matrices(const matrix &a_matrix, const matrix &b_matrix) {
    matrix s_matrix((a_matrix.size()), std::vector<double>(a_matrix.size()));
    for (size_t i = 0; i < a_matrix.size(); ++i)
        for (size_t j = 0; j < a_matrix[0].size(); ++j)
            s_matrix[i][j] = a_matrix[i][j] - b_matrix[i][j];
    return s_matrix;
}

/*
 * Multiplies matrix A by matrix B
 */
matrix multiply_matrices(const matrix &a_matrix, const matrix &b_matrix) {
    matrix m_matrix((a_matrix.size()), std::vector<double>(a_matrix.size()));
    for (size_t i = 0; i < a_matrix.size(); ++i) {
        for (size_t j = 0; j < a_matrix[0].size(); ++j) {
            for (size_t k = 0; k < a_matrix.size(); k++)
                m_matrix[i][j] += a_matrix[i][k] * b_matrix[k][j];
        }
    }
    return m_matrix;
}

/*
 * Multiplies matrix A by scalar value
 */
matrix multiply_by_scalar(const matrix &a_matrix, double val) {
    matrix m_matrix((a_matrix.size()), std::vector<double>(a_matrix.size()));
    for (size_t i = 0; i < a_matrix.size(); ++i)
        for (size_t j = 0; j < a_matrix[0].size(); ++j)
            m_matrix[i][j] = a_matrix[i][j] * val;
    return m_matrix;
}

/*
 * Calculate matrix trace - sum of elements on diagonal.
 */
double calculate_matrix_trace(const matrix &mtx) {
    double trace = 0;
    for (int i = 0; i < mtx.size(); ++i)
        for (int j = 0; j < mtx[0].size(); ++j)
            if (i == j)
                trace += mtx[i][j];
    return trace;
}

/*
 * Calculates residual matrix, which measures how matrix B
 * differs from matrix A.
 */
matrix calculate_R_matrix(const matrix &I_matrix, const matrix &BxA_matrix) {
    matrix R_matrix = substract_matrices(I_matrix, BxA_matrix);
    return R_matrix;
}

/*
 * If value of R_matrix is small enough, return true, else return false.
 */
bool check_R_matrix(const matrix &R_matrix) {
    double R_matrix_sum = 0.0;
    for (const auto &raw : R_matrix)
        for (const auto &element : raw)
            R_matrix_sum += fabs(element);

    return R_matrix_sum < EPSILON_RESIDUAL;
}

/*
 * Calculates first aproximation of inversed A_matrix called B_matrix
 */
matrix calculate_first_B_matrix(const matrix &trans_A_matrix,
                                const matrix &A_matrix) {
    matrix B_matrix = multiply_by_scalar(
        trans_A_matrix, 1 / calculate_matrix_trace(
                                multiply_matrices(trans_A_matrix, A_matrix)));
    return B_matrix;
}

/*
 * Calculates next aproximation of inversed A_matrix called B_matrix
 */
matrix calculate_next_B_matrix(const matrix &B_matrix,
                               const matrix &BxA_matrix) {
    auto f1 = std::async(std::launch::async, multiply_by_scalar, B_matrix, 2);
    auto f2 =
        std::async(std::launch::async, multiply_matrices, BxA_matrix, B_matrix);
    matrix next_B_matrix = substract_matrices(f1.get(), f2.get());
    return next_B_matrix;
}

/*
 * Inverses matrix with fast iterative method
 */
matrix inverse_matrix_iterative(const matrix &A_matrix) {
    std::future<matrix> f1 =
        std::async(std::launch::async, transpose_matrix, A_matrix);
    std::future<matrix> f2 =
        std::async(std::launch::async, create_identity_matrix, A_matrix.size());
    matrix B_matrix = calculate_first_B_matrix(f1.get(), A_matrix);
    const matrix I_matrix = f2.get();
    matrix R_matrix =
        calculate_R_matrix(I_matrix, multiply_matrices(B_matrix, A_matrix));

    while (!check_R_matrix(R_matrix)) {
        matrix BxA_matrix = multiply_matrices(B_matrix, A_matrix);
        std::future<matrix> f_B = std::async(
            std::launch::async, calculate_next_B_matrix, B_matrix, BxA_matrix);
        std::future<matrix> f_R = std::async(
            std::launch::async, calculate_R_matrix, I_matrix, BxA_matrix);
        R_matrix = f_R.get();
        B_matrix = f_B.get();
    }

    return B_matrix;
}

#endif
