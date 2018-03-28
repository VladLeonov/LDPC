#include "matrix.h"
#include "ldpc.h"

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
