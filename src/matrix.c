#include <stdlib.h>

#include "matrix.h"


#define TRUE !0
#define FALSE 0


int calculate_sum_row_elements(matrix G, int row_index) {
    int i = 0;
    int sum = 0;
    
    for (i = 0; i < G.columns; i++) {
        sum += G.body[row_index][i];
    }
    
    return sum;
}


float calculate_sum_coloumn_elements(int coloumns, float array[][coloumns],
                                     int coloumn_index, int rows) {
    int i = 0;
    float sum = 0.0;
    
    for (i = 0; i < rows; i++) {
        sum += array[i][coloumn_index];
    }
    
    return sum;
}


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


matrix create_unit_matrix(int rows) {
    matrix M = create_zero_matrix(rows, rows);
    int i = 0;
    
    for (i = 0; i < rows; i++) {
        M.body[i][i] = 1;
    }
    
    return M;
}

matrix create_void_matrix() {
    matrix M;
    M.body = NULL;
    M.rows = 0;
    M.columns = 0;
    
    return M;
}

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


int is_void_matrix(matrix M) {
    return !((M.rows > 0) && (M.columns > 0));
}


void free_matrix(matrix M) {
    int i = 0;
    
    for (i = 0; i < M.rows; i++) {
        free(M.body[i]);
    }
    free(M.body);
    
    M = create_void_matrix();
}


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
