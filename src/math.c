#include "matrix.h"
#include "ldpc.h"
#include "stdlib.h"
#include "math.h"

#define MAXITER 10000
int* gauss_elimination(matrix G) {
	int k = G.rows;
	int n = G.columns;
    int* information_set = (int*) malloc(k * sizeof(int));
    int i, j;
    for (i = 0; i < k; i++){
        information_set[i] = -1;
    }

    i = 0;
    int column_number = 0;

    while (i < k) {

        //throw away all-zero row
        while (sum_rows(G, i) == 0) {
            for (j = i; j < (k - 1); j++) {
                G.body[j] = G.body[j + 1];
            }
            k--;
            G.rows--;

            information_set[k] = -1;

            if (i >= k) {
                break;
            }
        }

        if (i >= k) {
            break;
        }

        //find 1 in the column
        int b = 0;
        j = i - 1;

        while ((b == 0) && (j < k - 1)) {
            j++;
            b = G.body[j][column_number];
        }

        if (b == 1) {
            information_set[i] = column_number;

               // remove other ones from this column
               int o, u;
               for (o = 0; o < k; o++) {
                   if (o == j) continue;

                   if (G.body[o][column_number]==1) {
                       for (u = 0; u < n; u++) {
                           G.body[o][u] ^= G.body[j][u];
                    }
                }
            }

            // replace rows i and j
            int *buff = G.body[i];
            G.body[i] = G.body[j];
            G.body[j] = buff;

            i++;
        }

        column_number++;  // keep the same row, another column
    }
    return information_set;
}

int sum_rows(matrix G, int row_index) {
    int i;
    int sum = 0;
    for (i = 0; i < G.columns; i++) {
        sum += G.body[row_index][i];
    }
    return sum;
}

void fill_with_permutation(int *x, int n) {
	int i, j, temp;
    for (i = 0; i < n; i++) {
        x[i] = i;
    }
    for (i = n-1; i > 0; i--) {
        j = rand() % (i+1);
        temp = x[j];
        x[j] = x[i];
        x[i] = temp;
    }
}

matrix create_G_from_H_matrix(matrix H, columns_metadata columns_mdata) {
	int i, j;
	matrix G = create_zero_matrix(H.columns - H.rows, H.columns);
    matrix P = create_empty_matrix(H.rows, columns_mdata.information_size);
    for (j = 0; j < columns_mdata.information_size; j++) {
        for (i = 0; i < H.rows; i++) {
            P.body[i][j] = H.body[i][columns_mdata.information_set[j]];
        }
    }

    matrix PT = transpose_matrix(P);
    for (j = 0; j < columns_mdata.check_size; j++) {
        for (i = 0; i < G.rows; i++) {
            G.body[i][columns_mdata.check_set[j]] = PT.body[i][j];
        }
    }

    matrix U = create_unit_matrix(G.rows);
    for (j = 0; j < columns_mdata.information_size; j++) {
        for (i = 0; i < G.rows; i++) {
            G.body[i][columns_mdata.information_set[j]] = U.body[i][j];
        }
    }

    free_matrix(P);
    free_matrix(PT);
    free_matrix(U);

    return G;
}

int get_indexes_of_common_elements(int *arr_a, int *arr_b, int *result, int len_a, int len_b) {
    int i, j, result_length;
    j = 0;
    result_length = 0;

    for (i = 0; i < len_a; i++) {
        for (j = 0; j < len_b; j++) {
            if (arr_a[i] == arr_b[j]) {
                result[result_length] = i;
                result_length++;
                break;
            } else if (arr_a[i] < arr_b[j]) {
                break;
            } else {
                continue;
            }
        }
    }

    return result_length;
}

matrix  get_hard_from_soft(double soft[], int length) {
    matrix result = create_zero_matrix(1, length);
    
    int i;
    for (i = 0; i < length; i++) {
        if (soft[i] > 0) {
            result.body[0][i] = 1;
        } else {
            result.body[0][i] = 0;
        }
    }

    return result;
}

double log_tahn(double value) {
    double t = exp(value);
    return log((t + 1)/(t - 1));
}

double max(double value1, double value2) {
    if (value1 > value2) {
        return value1;
    } else {
        return value2;
    }
}

double min(double value1, double value2) {
    if (value1 > value2) {
        return value2;
    } else {
        return value1;
    }
}

double sum_array(double *array, int length) {
    int i;
    double result = 0.0;
    for (int i = 0; i < length; i++) {
        result += array[i];
    }
    return result;
}

int decode_belief_propogandation(ldpc ldpc_object, double *y, matrix *hard_solution) {
    matrix H = ldpc_object.H;
    int r = H.rows;
    int n = H.columns;
    double Z[r][n];
    double L[r][n];
    non_zero_data C = ldpc_object.C;
    non_zero_data V = ldpc_object.V;

    int i, j;
    for (i = 0; i < r; i++) {
        for (j = 0; j < n; j++) {
            Z[i][j] = 0;
            L[i][j] = 0;
        }
    }

    double soft[n];
    for (i = 0; i < n; i++) {
        soft[i] = y[i];
    }

    matrix hard = get_hard_from_soft(soft, n);

    matrix syndrome = multiply_matrices(hard, transpose_matrix(H));
    int iter = 0;
    int L_element = 0;
    int index_C = 0;
    int k;
    
    while ((iter < MAXITER) && (sum_syndrome(syndrome) != 0)) {
        //H columns processing
        double sum_Z;
        for (i = 0; i < n; i++) {
        	sum_Z = 0;
            for (j = 0; j < C.element_length[i]; j++) {
                sum_Z += Z[C.element_data[i][j]][i];
            }
            for (j = 0; j < C.element_length[i]; j++) {
                index_C = C.element_data[i][j];
                L[index_C][i] = y[i] + sum_Z - Z[index_C][i];
            }
        }
        
        //H rows processing
        double a, b;
        for (i = 0; i < n; i++) {
            for (j = 0; j < C.element_length[i]; j++) {
                int index_C = C.element_data[i][j];
                a = 1;
                b = 0;
                for (k = 0; k < V.element_length[index_C]; k++) {
                    L_element = L[index_C][V.element_data[index_C][k]];
                    a *= sign(L_element);
                    b += log_tahn(L_element);
                }
                a *= sign(L[index_C][i]);
                b -= log_tahn(L[index_C][i]);
                Z[index_C][i] = max(min(a * log_tahn(b), 19.07), -19.07);
            }
        }

        //result forming
        for (i = 0; i < (n - r); i++) {
            soft[i] = y[i] + sum_array(Z[i], n);
        }
        hard = get_hard_from_soft(soft, n);
        free_matrix(syndrome);
        syndrome = multiply_matrices(hard, transpose_matrix(H));
        iter++;

    }
    
    hard_solution->body = hard.body;
    hard_solution->columns = hard.columns;
    hard_solution->rows = hard.rows;
    
    free_matrix(syndrome);
    
    return iter;
}

double* product_matrix(double **data, int rows, int columns) {
    int i, j;
    double *result = (double*)malloc(columns * sizeof(double));

    for (i = 0; i < columns; i++) {
        result[i] = 1;
    }

    for (i = 0; i < columns; i++) {
        for (j = 0; j < rows; j++) {
            result[i] *= data[i][j];
        }
    }

    return result;
}

double* sum_matrix_columns(int rows, int columns, double data[rows][columns]) {
    int i, j;
    double result[columns];

    for (i = 0; i < columns; i++) {
        result[i] = 0;
    }

    for (i = 0; i < columns; i++) {
        for (j = 0; j < rows; j++) {
            result[i] += data[j][i];
        }
    }

    return result;
}

int sum_syndrome(matrix syndrome) {
    int result = 0;
    int i = 0;
    for (int i = 0; i < syndrome.columns; i++) {
        result += syndrome.body[0][i];
    }

    return result;
}

int sign(double value) {
    if (value > 0) {
        return 1;
    } else if (value < 0) {
        return -1;
    } else {
        return 0;
    }
}
