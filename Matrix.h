#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <cmath>
#include <iostream>
#include <random>
#include <string>
#include <vector>
#include <mpi.h>

using matrix = std::vector<std::vector<double>>;
const double EPSILON = 1e-11;
const int DIMENSION = 1000;
const std::pair<int, int> RANDOM_RANGE = {1, 10000};

class ZeroDiagonalEception : std::exception {};

// global mpi array
double *mpi_array;

void print_matrix(const matrix &mtx);
matrix create_random_matrix(int dim);
matrix create_identity_matrix(int dim);
matrix transpose_matrix(const matrix &mtx);
matrix substract_matrices(const matrix &a_matrix, const matrix &b_matrix);
void multiply_matrices_partially(const matrix &a_matrix, const matrix &b_matrix,
                                 int start, int end);
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
void calculate_column(int k, matrix &A, matrix &X);
void decompose_to_LU(matrix &a_matrix);
matrix inverse_matrix(matrix &mtx);

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
    for (auto &row : random_matrix)
        for (auto &element : row)
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
 * Multiplies part of matrix A by part of matrix B
 */
void multiply_matrices_partially(const matrix &a_matrix, const matrix &b_matrix,
                                 int start, int end) {
    size_t size = a_matrix.size();
    mpi_array = new double [(end-start) * size];
    int n = 0;
    for (size_t i = start; i < end; ++i) {
        for (size_t j = 0; j < size; ++j) {
            for (size_t k = 0; k < size; k++)
                mpi_array[n] += a_matrix[i][k] * b_matrix[k][j];
	    ++n;
        }
    }
}

/*
 * Multiplies matrix A by matrix B - parallel
 */
matrix multiply_matrices(const matrix &a_matrix, const matrix &b_matrix) {
    MPI_Barrier(MPI_COMM_WORLD);
    size_t m_size = a_matrix.size();
    int proc_index, proc_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_index);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_size);
    MPI_Status status;
    int rest = m_size % proc_size;
    int quotient = std::floor(m_size / proc_size);
    int start = proc_index * quotient;
    int end = (proc_index + 1) * quotient;
    std::cout << "start: " << start << " end: " << end << std::endl;

    matrix m_matrix((m_size), std::vector<double>(m_size));
    multiply_matrices_partially(a_matrix, b_matrix, start, end);

    if (proc_index != 0)
	MPI_Ssend(mpi_array, quotient * m_size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);

    if (proc_index == 0) {
	int n = 0;
	for (size_t i = start; i < end; ++i) {
	    for (size_t j = 0; j < m_size; ++j) {
		m_matrix[i][j] = mpi_array[n];
	    }
	    ++n;
	}
	for (int source = 1; source < proc_size; ++source) {
	    MPI_Recv(mpi_array, quotient * m_size, MPI_DOUBLE, source, 0,
		    MPI_COMM_WORLD, &status);
	    int n = 0;
	    for (size_t i = start; i < end; ++i) {
		for (size_t j = 0; j < m_size; ++j) {
		    m_matrix[i][j] = mpi_array[n];
		}
		++n;
	    }
	}
    }

    // if (rest != 0)
	// multiply_matrices_partially(a_matrix, b_matrix, m_matrix,
				    // proc_size * std::floor(m_size / proc_size),
				    // m_size);
    std::cout << "dochodze tuu: " << proc_index << std::endl;
    delete mpi_array;
    mpi_array = nullptr;
    return m_matrix;
}

/*
 * Multiplies matrix A by matrix B - no parallelisation
 */
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
    // double R_matrix_sum = 0.0;
    for (const auto &row : R_matrix)
        for (const auto &element : row)
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
    matrix next_B_matrix = substract_matrices(multiply_by_scalar(B_matrix, 2),
					      multiply_matrices(BxA_matrix,
								B_matrix));
    return next_B_matrix;
}

/*
 * Inverses matrix with fast iterative method
 */
matrix inverse_matrix_iterative(const matrix &A_matrix) {
    int proc_index, proc_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_index);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_size);
    auto T_matrix = transpose_matrix(A_matrix);

    matrix B_matrix = calculate_first_B_matrix(T_matrix, A_matrix);
    auto I_matrix = create_identity_matrix(A_matrix.size());
    matrix R_matrix =
        calculate_R_matrix(I_matrix, multiply_matrices(B_matrix, A_matrix));

    MPI_Barrier(MPI_COMM_WORLD);
    while (true) {
        matrix BxA_matrix = multiply_matrices(B_matrix, A_matrix);
	MPI_Barrier(MPI_COMM_WORLD);
	if (proc_index == 0)
	    R_matrix = calculate_R_matrix(I_matrix, BxA_matrix);
	MPI_Barrier(MPI_COMM_WORLD);

        B_matrix = calculate_next_B_matrix(B_matrix, BxA_matrix);
	MPI_Barrier(MPI_COMM_WORLD);
	if (!check_R_matrix(R_matrix)) {

	}

    }
    return B_matrix;
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


/*** Matrix inversion using non-iterative method. ***/
/*** Used for tests purposes. ***/

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
#endif
