/**
    LDPC
    matrix.c
    Purpose: Ñontains subsidary operations with matrices
	for using by other methods.

    @author Leonov V.R.
    @version 08.05.18
*/


#include <stdlib.h>

#include "matrix.h"


#define TRUE !0
#define FALSE 0


//Calculates the sum of all elements of a matrix row.
int calculate_sum_row_elements(matrix G, int row_index) {
    int i = 0;
    int sum = 0;
    
    for (i = 0; i < G.columns; i++) {
        sum += G.body[row_index][i];
    }
    
    return sum;
}


//Calculates the sum of all elements of a two-dimensional array coloumn.
float calculate_sum_column_elements(int columns, int rows, 
                                    float array[rows][columns],
									int column_index) {
    int i = 0;
    float sum = 0.0;
    
    for (i = 0; i < rows; i++) {
        sum += array[i][column_index];
    }
    
    return sum;
}


//Creates matrix that randomly filled with zeros and ones.
matrix create_random_matrix(int rows, int columns) {
    int i = 0, j = 0;
    matrix M = create_empty_matrix(rows, columns);
    
    for (i = 0; i < rows; i++) {
    	for (j = 0; j < columns; j++) {
    		M.body[i][j] = rand() % 2;
    	}
    }
    
    return M;
}


//Creates matrix with uninitialized elements.
matrix create_empty_matrix(int rows, int columns) {
    matrix M;
    M.rows = rows;
    M.columns = columns;
    M.body = (int**) malloc(rows * sizeof(int*));
    
    int i = 0;
    
    for (i = 0; i < rows; i++) {
        M.body[i] = (int*) malloc(columns * sizeof(int));
    }
    
    return M;
}


//Creates zero matrix.
matrix create_zero_matrix(int rows, int columns) {
    matrix M;
    M.rows = rows;
    M.columns = columns;
    M.body = (int**) malloc(rows * sizeof(int*));
    
    int i = 0;
    
    for (i = 0; i < rows; i++) {
        M.body[i] = (int*) calloc(columns, sizeof(int));
    }
    
    return M;
}


//Creates unit matrix.
matrix create_unit_matrix(int rows) {
    matrix M = create_zero_matrix(rows, rows);
    int i = 0;
    
    for (i = 0; i < rows; i++) {
        M.body[i][i] = 1;
    }
    
    return M;
}


//Creates matrix without elements.
matrix create_void_matrix() {
    matrix M;
    M.body = NULL;
    M.rows = 0;
    M.columns = 0;
    
    return M;
}


//Multiplies two matrices.
matrix multiply_matrices(matrix M1, matrix M2) {
    if (M1.columns != M2.rows) return create_void_matrix();
    
    matrix M = create_zero_matrix(M1.rows, M2.columns);
    int i = 0, j = 0, k = 0;
    
    for (i = 0; i < M1.rows; i++) {
        for (j = 0; j < M2.columns; j++) {
            for (k = 0; k < M1.columns; k++) {
                M.body[i][j] ^= M1.body[i][k] * M2.body[k][j];
            }
        }
    }
    
    return M;
}


//Transposes matrix.
matrix transpose_matrix(matrix M) {
    matrix MT = create_empty_matrix(M.columns, M.rows);
    int i = 0, j = 0;
    
    for (i = 0; i < M.rows; i++) {
        for (j = 0; j < M.columns; j++) {
            MT.body[j][i] = M.body[i][j];
        }
    }
    
    return MT;
}


//Ñhecks whether the matrix contains elements.
int is_void_matrix(matrix M) {
    return !((M.rows > 0) && (M.columns > 0));
}


//Frees the memory occupied by the matrix, makes it an void matrix.
void free_matrix(matrix M) {
    int i = 0;
    
    for (i = 0; i < M.rows; i++) {
        free(M.body[i]);
    }
    free(M.body);
    
    M = create_void_matrix();
}


//Makes copy of matrix.
matrix copy_matrix(matrix M){
    matrix copy_M = create_empty_matrix(M.rows, M.columns);
    int i = 0, j = 0;
    
    for (i = 0; i < M.rows; i++) {
        for (j = 0; j < M.columns; j++) {
            copy_M.body[i][j] = M.body[i][j];
        }
    }
    
    return copy_M;
}


//Makes copy of matrix part.
matrix copy_matrix_part(matrix M, int rows, int columns) {
	int num_rows_to_copy = 0, num_columns_to_copy = 0;

	if (M.rows > rows) {
		num_rows_to_copy = rows;
	} else {
		num_rows_to_copy = M.rows;
	}

	if (M.columns > columns) {
        num_columns_to_copy = columns;
	} else {
        num_columns_to_copy = M.columns;
	}

	matrix part_M = create_zero_matrix(rows, columns);
	int i = 0, j = 0;

	for (i = 0; i < num_rows_to_copy; i++){
		for (j = 0; j < num_columns_to_copy; j++) {
			part_M.body[i][j] = M.body[i][j];
		}
	}

	return part_M;
}


//Compares two matrices.
char compare_matrices(matrix M1, matrix M2) {
    if ((M1.columns != M2.columns) || (M1.rows != M2.rows)) return FALSE;

    int i = 0, j = 0;
    
    for (i = 0; i < M1.rows; i++) {
        for (j = 0; j < M1.columns; j++) {
            if (M1.body[i][j] != M2.body[i][j]) return FALSE;
        }
    }

    return TRUE;
}


//Converts two-dimensional array to matrix.
matrix array_to_matrix(int rows, int columns, int array[rows][columns]) {
    matrix new_matrix = create_empty_matrix(rows, columns);
    int i, j;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < columns; j++) {
            new_matrix.body[i][j] = array[i][j];
        }
    }
    return new_matrix;
}
