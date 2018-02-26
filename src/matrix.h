#ifndef MATRIX
#define MATRIX

typedef struct {
	char **body;
	int rows, columns;
} matrix;

matrix create_empty_matrix(int rows, int columns);
matrix create_zero_matrix(int rows, int columns);
matrix create_unit_matrix(int rows);
matrix create_sparse_matrix(int rows, int columns);

matrix multiply_matrices(matrix M1, matrix M2);
matrix transpose_matrix(matrix M);
matrix combine_matrices(matrix M1, matrix M2);

matrix create_void_matrix();
int is_void_matrix(matrix M);

void free_matrix(matrix M);

#endif
