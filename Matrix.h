#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <cmath>
#include <future>
#include <iostream>
#include <random>
#include <string>
#include <vector>


using matrix = std::vector<std::vector<double>>;
const double EPSILON = 1e-11;
const int DIMENSION = 1000;
const std::pair<int, int> RANDOM_RANGE = {1, 10000};

class ZeroDiagonalEception : std::exception {};

matrix create_identity_matrix(int dim);
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
void multiply_matrices_partially(const matrix &a_matrix, const matrix &b_matrix, matrix &output, int start, int end);

/*
 * Prints out the matrix.
 * Use to debug.
 */
void print_matrix(const matrix &mtx) {
    for (const auto &row : mtx) {
        for (const auto &element : row)
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
 * Transposes square matrix.
 */
matrix transpose_matrix(const matrix &mtx) {
    size_t size = mtx.size();
    matrix t_matrix((size), std::vector<double>(size));
    for (size_t i = 0; i < size; ++i)
        for (size_t j = 0; j < size; ++j)
            t_matrix[j][i] = mtx[i][j];
    return t_matrix;
}

/*
 * Substracts matrix A from matrix B
 */
matrix substract_matrices(const matrix &a_matrix, const matrix &b_matrix) {
    size_t size = a_matrix.size();
    matrix s_matrix((size), std::vector<double>(size));
    for (size_t i = 0; i < size; ++i)
        for (size_t j = 0; j < size; ++j)
            s_matrix[i][j] = a_matrix[i][j] - b_matrix[i][j];
    return s_matrix;
}

/*
 * Multiplies matrix A by matrix B
 */
void multiply_matrices_partially(const matrix &a_matrix, const matrix &b_matrix, matrix &m_matrix, int start, int end){
  size_t size = a_matrix.size();
  for (size_t i = start; i < end; ++i) {
      for (size_t j = 0; j < size; ++j) {
          for (size_t k = 0; k < size; k++)
              m_matrix[i][j] += a_matrix[i][k] * b_matrix[k][j];
      }
  }
}

matrix multiply_matrices(const matrix &a_matrix, const matrix &b_matrix) {
    size_t m_size = a_matrix.size();
    matrix m_matrix((m_size), std::vector<double>(m_size));

    size_t thread_size = 10;

    std::vector<std::future<void>> v;

    for (auto i=0; i < thread_size; ++i){
      v.push_back(std::async(std::launch::async, multiply_matrices_partially, std::ref(a_matrix), std::ref(b_matrix), std::ref(m_matrix), i * m_size/thread_size, m_size/thread_size * (i+1)));
    }

    for (auto it = v.begin(); v.end() != it; it++) {
        (*it).get();
      }

    return m_matrix;
}


// matrix multiply_matrices(const matrix &a_matrix, const matrix &b_matrix) {
//     matrix m_matrix((a_matrix.size()), std::vector<double>(a_matrix.size()));
//     for (size_t i = 0; i < a_matrix.size(); ++i) {
//         for (size_t j = 0; j < a_matrix[0].size(); ++j) {
//             for (size_t k = 0; k < a_matrix.size(); k++)
//                 m_matrix[i][j] += a_matrix[i][k] * b_matrix[k][j];
//         }
//     }
//     return m_matrix;
// }

/*
 * Multiplies matrix A by scalar value
 */
matrix multiply_by_scalar(const matrix &a_matrix, double val) {
    size_t size = a_matrix.size();
    matrix m_matrix((size), std::vector<double>(size));

    for (size_t i = 0; i < size; ++i)
        for (size_t j = 0; j < size; ++j)
            m_matrix[i][j] = a_matrix[i][j] * val;
    return m_matrix;
}

/*
 * Calculate matrix trace - sum of elements on diagonal.
 */
double calculate_matrix_trace(const matrix &mtx) {
    double trace = 0;
    size_t size = mtx.size();
    for (int i = 0; i < size; ++i)
        for (int j = 0; j < size; ++j)
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
    //double R_matrix_sum = 0.0;
    for (const auto &raw : R_matrix)
        for (const auto &element : raw)
            if (fabs(element) > EPSILON)
              return false;

    return true;
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
