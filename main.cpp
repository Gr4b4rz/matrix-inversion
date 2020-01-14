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

    MPI_Init(&argc, &argv);

    matrix random_matrix = create_random_matrix(dim);
    matrix m = inverse_matrix_iterative(random_matrix);
    int proc_index, proc_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_index);

    if (proc_index == 0)
	std::cout << "Matrix has been successfully inversed" << std::endl;

    MPI_Finalize();

    return 0;
}

