/**
    LDPC
    matrix.h
    Purpose: Ñontains subsidary operations with matrices
	for using by other methods.

    @author Leonov V.R.
    @version 24.04.18
*/


#ifndef MATRIX_CORRECT_VERSION_H_INCLUDED
#define MATRIX_CORRECT_VERSION_H_INCLUDED


//Matrix of integers.
typedef struct {
    int **body;		//Double array with matrix elements.
    int rows;		//Number of rows in the matrix.
	int columns;	//Number of columns in the matrix.
} matrix;

/**
    Calculates the sum of all elements of a matrix row.

    @param G Matrix.
    @param row_index Index of the required row.
    @return Sum of elements of row.
*/
int calculate_sum_row_elements(matrix G, int row_index);

/**
    Calculates the sum of all elements of a two-dimensional array coloumn.

    @param coloumns Number of coloumns in array.
    @param array Two-dimensional array.
    @param coloumn_index Index of the required coloumn.
    @param rows Number of rows in array.
    @return Sum of elements of coloumn.
*/
float calculate_sum_coloumn_elements(int coloumns, float array[][coloumns], 
                                     int coloumn_index, int rows);

/**
    Creates matrix that randomly filled with zeros and ones.

    @param rows Number of rows in matrix.
    @param coloumns Number of coloumns in matrix.
    @return Generated matrix.
*/
matrix create_random_matrix(int rows, int columns);

/**
    Creates matrix with uninitialized elements.

    @param rows Number of rows in matrix.
    @param coloumns Number of coloumns in matrix.
    @return Generated matrix.
*/
matrix create_empty_matrix(int rows, int columns);

/**
    Creates zero matrix.

    @param rows Number of rows in matrix.
    @param coloumns Number of coloumns in matrix.
    @return Generated matrix.
*/
matrix create_zero_matrix(int rows, int columns);

/**
    Creates unit matrix.

    @param rows Number of rows/coloumns in matrix.
    @return Generated matrix.
*/
matrix create_unit_matrix(int rows);

/**
    Multiplies two matrices.

    @param M1 First matrix.
    @param M2 Second matrix.
    @return M1 * M2.
*/
matrix multiply_matrices(matrix M1, matrix M2);

/**
    Transposes matrix.

    @param M Matrix.
    @return Transposed M.
*/
matrix transpose_matrix(matrix M);

/**
    Ñhecks whether the matrix contains elements.

    @param M Matrix.
    @return FALSE, if matrix sizes not zero, TRUE otherwise.
*/
int is_void_matrix(matrix M);

/**
    Frees the memory occupied by the matrix, makes it an void matrix.

    @param M Matrix.
*/
void free_matrix(matrix M);

/**
    Makes copy of matrix.

    @param M Matrix.
    @return Matrix, that equal M.
*/
matrix copy_matrix(matrix M);

/**
    Makes copy of matrix part.

    @param M Matrix.
    @param rows Number of rows in matrix part.
    @param coloumns Number of coloumns in matrix part.
    @return Matrix, that contain part of M.
*/
matrix copy_matrix_part(matrix M, int rows, int columns);

/**
    Compares two matrices.

    @param M1 First matrix.
    @param M2 Second matrix.
    @return TRUE, if matrices equal, FALSE otherwise.
*/
char compare_matrices(matrix M1, matrix M2);

#endif // MATRIX_CORRECT_VERSION_H_INCLUDED
