#include "Matrix.h"

int main(int argc, char **argv) {
    int dim;
    if (argc != 2 || atoi(argv[1]) < 1)
        dim = DIMENSION;
    else
        dim = atoi(argv[1]);

    matrix random_matrix = create_random_matrix(dim);
    inverse_matrix_iterative(random_matrix);
    // inverse_matrix(random_matrix);
    std::cout << "Matrix has been successfully inversed" << std::endl;

    return 0;
}
