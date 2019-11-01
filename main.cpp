#include <cmath>
#include <iostream>
#include <random>
#include <string>
#include <vector>

using matrix = std::vector<std::vector<double>>;
const double EPSILON = 1e-16;
const int DIMENSION = 1000;
const std::pair<int, int> RANDOM_RANGE = {1, 10000};

class ZeroDiagonalEception : std::exception {};

void calculate_column(int k, matrix &A, matrix &X);
void decompose_to_LU(matrix &a_matrix);
matrix create_identity_matrix();
matrix inverse_matrix(matrix &mtx);
void print_matrix(const matrix &mtx);
matrix create_random_matrix(int dim);


int main(int argc, char **argv) {
    int dim;
    if (argc != 2 || atoi(argv[1]) < 1)
        dim = DIMENSION;
    else
        dim = atoi(argv[1]);

    matrix random_matrix = create_random_matrix(dim);
    // print_matrix(random_matrix);
    auto reverse_matrix = inverse_matrix(random_matrix);
    // print_matrix(reverse_matrix);
    std::cout << "Matrix has been successfully inversed" << std::endl;

    return 0;
}

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
matrix create_identity_matrix(const int dim) {
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
