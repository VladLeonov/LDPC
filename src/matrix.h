#ifndef MATRIX_CORRECT_VERSION_H_INCLUDED
#define MATRIX_CORRECT_VERSION_H_INCLUDED

typedef struct {
    int **body;
    int rows, columns;
} matrix;

int calculate_sum_row_elements(matrix G, int row_index);
float calculate_sum_coloumn_elements(int coloumns, float array[][coloumns], int coloumn_index, int rows);
matrix create_random_matrix(int rows, int columns);
matrix create_empty_matrix(int rows, int columns);
matrix create_zero_matrix(int rows, int columns);
matrix create_unit_matrix(int rows);
matrix multiply_matrices(matrix M1, matrix M2);
matrix transpose_matrix(matrix M);
matrix create_void_matrix();
int is_void_matrix(matrix M);
void free_matrix(matrix M);
matrix copy_matrix(matrix matrix_object);
matrix copy_matrix_part(matrix old_matrix_object, int rows, int columns);
char compare_matrices(matrix matrix1, matrix matrix2);

#endif // MATRIX_CORRECT_VERSION_H_INCLUDED
