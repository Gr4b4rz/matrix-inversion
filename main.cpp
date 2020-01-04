#include <future>
#include <iostream>
#include <string>
#include <thread>
#include <mpi.h>
#include "Matrix.h"

void solicit(int *pN, int *pnprocs, int mynum);

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
    int mynum, nrprocs, N = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &mynum);
    MPI_Comm_size(MPI_COMM_WORLD, &nrprocs);
    // size is always 1?; Weird ...; Check if u got the same
    std::cout << "Comm_rank=" << mynum << ",Comm_size=" << nrprocs << std::endl;
    // solicit(&N, &nrprocs, mynum);
    // std::cout << "N: " << N << std::endl;
    MPI_Finalize();

    return 0;
}

/*
 * Get a value for N, the number of intervals in the approximation.
 * (Parallel versions: master instance reads in N and then
 * broadcasts N to all the other instances of the program.)
 */
void solicit (int *pN, int *pnprocs, int mynum)
{
    int source = 0;
    if (mynum == 0) {
	printf ("Enter number of approximation intervals:(0 to exit)\n");
	scanf ("%d", pN);
    }
    MPI_Bcast(pN, 1, MPI_INT, source, MPI_COMM_WORLD);
}

