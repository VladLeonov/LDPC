#include "matrix.h"

#include <malloc.h>
#include <time.h>
#include <stdlib.h>

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
    /*matrix M = create_empty_matrix(rows, columns);
    int i, j;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < columns; j++) {
            M.body[i][j] = 0;
        }
    }*/
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


matrix create_sparse_matrix(int rows, int columns) {
    srand(time(NULL));
    matrix M = create_zero_matrix(rows, columns);
    int i;
    for (i = 0; i < rows; i++) {
        M.body[i][rand() % columns] = 1;
    }
    for (i = 0; i < columns; i++) {
        M.body[rand() % rows][i] = 1;
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


matrix combine_matrices(matrix M1, matrix M2) {
    if (M1.rows != M2.rows) {
        return create_void_matrix();
    }
    matrix M = create_empty_matrix(M1.rows, M1.columns + M2.columns);
    int i, j;
    for (i = 0; i < M1.rows; i++) {
        for (j = 0; j < M1.columns; j++) {
            M.body[i][j] = M1.body[i][j];
        }
        for (j = 0; j < M2.columns; j++) {
            M.body[i][j + M1.columns] = M2.body[i][j];
        }
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


int is_void_matrix(matrix M) {
    return ((M.rows > 0) && (M.columns > 0));
}


void free_matrix(matrix M) {
    int i;
    for (i = 0; i < M.rows; i++) {
        free(M.body[i]);
    }
    free(M.body);
    M = create_void_matrix();
}

void print_matrix(matrix M) {
    int i, j;
    for (i = 0; i < M.rows; i++) {
        for (j = 0; j < M.columns; j++) {
            printf("%d ", M.body[i][j]);
        }
        printf("\n");
    }
}

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

non_zero_data get_non_zero_column_data(matrix matrix_object) {
    non_zero_data result;
    result.element_data = (int **)malloc(matrix_object.columns * sizeof(int *));
    result.element_length = (int *)malloc(matrix_object.columns * sizeof(int));

    int i, j, num_ones = 0;
    int buf_column[matrix_object.rows];
    for (i = 0; i < matrix_object.columns; i++) {
        for (j = 0; j < matrix_object.rows; j++) {
            if (matrix_object.body[j][i] != 0) {
                buf_column[num_ones] = j;
                num_ones ++;
            }
        }
        result.element_length[i] = num_ones;
        result.element_data[i] = (int *)malloc(num_ones * sizeof(int));

        for (j = 0; j < num_ones; j++) {
            result.element_data[i][j] = buf_column[j];
        }
        num_ones = 0;
    }

    return result;
}
