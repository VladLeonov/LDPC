#include "matrix.h"
#include "stdlib.h"

#define TRUE !0
#define FALSE 0

int sum_row_elements(matrix G, int row_index) {
    int i;
    int sum = 0;
    for (i = 0; i < G.columns; i++) {
        sum += G.body[row_index][i];
    }
    return sum;
}

float sum_coloumn_elements(int coloumns, float array[][coloumns], int coloumn_index, int rows) {
    int i;
    float result = 0.0;
    for (i = 0; i < rows; i++) {
        result += array[i][coloumn_index];
    }
    return result;
}

matrix create_random_matrix(int rows, int columns) {
    int i, j;
    matrix result = create_empty_matrix(rows, columns);
    for (i = 0; i < rows; i++) {
    	for (j = 0; j < columns; j++) {
    		result.body[i][j] = rand() % 2;
    	}
    }
    return result;
}

matrix create_empty_matrix(int rows, int columns) {
    matrix M;
    M.rows = rows;
    M.columns = columns;
    M.body = (int**) malloc(rows * sizeof(int*));
    int i;
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
    int i;
    for (i = 0; i < rows; i++) {
        M.body[i] = (int*) calloc(columns, sizeof(int));
    }
    return M;
}


matrix create_unit_matrix(int rows) {
    matrix M = create_zero_matrix(rows, rows);
    int i;
    for (i = 0; i < rows; i++) {
        M.body[i][i] = 1;
    }
    return M;
}


matrix multiply_matrices(matrix M1, matrix M2) {
    if (M1.columns != M2.rows) {
        return create_void_matrix();
    }
    matrix M = create_zero_matrix(M1.rows, M2.columns);
    int i, j, k;
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
    int i, j;
    for (i = 0; i < M.rows; i++) {
        for (j = 0; j < M.columns; j++) {
            MT.body[j][i] = M.body[i][j];
        }
    }
    return MT;
}


matrix create_void_matrix() {
    matrix M;
    M.body = NULL;
    M.rows = 0;
    M.columns = 0;
    return M;
}


int is_void_matrix(matrix M) {
    return !((M.rows > 0) && (M.columns > 0));
}


void free_matrix(matrix M) {
    int i;
    for (i = 0; i < M.rows; i++) {
        free(M.body[i]);
    }
    free(M.body);
    M = create_void_matrix();
}


matrix copy_matrix(matrix matrix_object){
    matrix new_matrix = create_empty_matrix(matrix_object.rows, matrix_object.columns);
    int i = 0, j = 0;
    for (i = 0; i < matrix_object.rows; i++) {
        for (j = 0; j < matrix_object.columns; j++) {
            new_matrix.body[i][j] = matrix_object.body[i][j];
        }
    }
    return new_matrix;
}

matrix copy_matrix_part(matrix old_matrix_object, int rows, int columns) {

	int i, j, new_rows, new_columns;

	if (old_matrix_object.rows > rows) {
		new_rows = rows;
	} else {
		new_rows = old_matrix_object.rows;
	}

	if (old_matrix_object.columns > columns) {
        new_columns = columns;
	} else {
        new_columns = old_matrix_object.columns;
	}

	matrix new_matrix_object = create_zero_matrix(rows, columns);

	for (i = 0; i < new_rows; i++){
		for (j = 0; j < new_columns; j++) {
			new_matrix_object.body[i][j] = old_matrix_object.body[i][j];
		}
	}

	return new_matrix_object;
}

char compare_matrices(matrix matrix1, matrix matrix2) {
    if ((matrix1.columns != matrix2.columns) || (matrix1.rows != matrix2.rows)) {
        return FALSE;
    }

    int i, j;
    for (i = 0; i < matrix1.rows; i++) {
        for (j = 0; j < matrix1.columns; j++) {
            if (matrix1.body[i][j] != matrix2.body[i][j]) {
                return FALSE;
            }
        }
    }

    return TRUE;
}
