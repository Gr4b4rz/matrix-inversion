#include <future>
#include <iostream>
#include <string>
#include <thread>
#include "Matrix.h"

void multiply_matrices1(const matrix &a_matrix, const matrix &b_matrix,
                        matrix &m_matrix, int start, int end);
matrix multiply_parallel(const matrix &a_matrix, const matrix &b_matrix);

int main(int argc, char **argv) {
    int dim;
    if (argc != 2 || atoi(argv[1]) < 1)
        dim = DIMENSION;
    else
        dim = atoi(argv[1]);

    std::cout << "Suggested number of threads = " <<
    std::thread::hardware_concurrency() << std::endl;
    matrix random_matrix = create_random_matrix(dim);
    matrix m = inverse_matrix_iterative(random_matrix);
    std::cout << "Matrix has been successfully inversed" << std::endl;

    return 0;
}
