#include <future>
#include <iostream>
#include <string>
#include <thread>
#include <mpi.h>
#include "Matrix.h"

int main(int argc, char **argv) {
    int dim;
    if (argc != 2 || atoi(argv[1]) < 1)
	dim = DIMENSION;
    else
	dim = atoi(argv[1]);

    matrix random_matrix = create_random_matrix(dim);
    matrix m = inverse_matrix_iterative(random_matrix);
    std::cout << "Matrix has been successfully inversed" << std::endl;

    MPI_Init(&argc, &argv);
    // MPI stuff
    MPI_Finalize();

    return 0;
}
