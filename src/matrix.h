#ifndef MATRIX
#define MATRIX

typedef struct {
    int **body;
    int rows, columns;
} matrix;

typedef struct{
    int **element_data;
    int *element_length;
}non_zero_data;

matrix create_empty_matrix(int rows, int columns);
matrix create_zero_matrix(int rows, int columns);
matrix create_unit_matrix(int rows);
matrix create_sparse_matrix(int rows, int columns);

matrix multiply_matrices(matrix M1, matrix M2);
matrix transpose_matrix(matrix M);
matrix combine_matrices(matrix M1, matrix M2);
matrix array_to_matrix(int rows, int columns, int array[rows][columns]);

matrix create_void_matrix();
int is_void_matrix(matrix M);

void free_matrix(matrix M);
void print_matrix(matrix M);
matrix copy_matrix(matrix matrix_object);
matrix copy_matrix_part(matrix old_matrix_object, int rows, int columns);

non_zero_data get_non_zero_column_data(matrix matrix_object);

char compare_matrices(matrix matrix1, matrix matrix2);
#endif
