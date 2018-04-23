#include <stdlib.h>
#include "ldpc_generator.h"
#include "matrix.h"
#include "subsidary_math.h"

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
        while (sum_row_elements(G, i) == 0) {
            for (j = i; j < (k - 1); j++) {
                G.body[j] = G.body[j + 1];
            }
            k--;

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
    for (i = k; i < G.rows; i++) {
    	G.body[i] = calloc(G.columns, sizeof(int));
	}
    return information_set;
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

ldpc create_ldpc(code_type type, int J, int K, int M) {

    matrix H = create_H_rand(type, J, K, M);
    matrix H_copy = copy_matrix(H);
    int* check_set = gauss_elimination(H_copy);
    columns_metadata columns_mdata = create_columns_metadata(check_set, K * M, J * M);

    matrix cutted_H;
    switch (type) {
        case Gallager:
            cutted_H = copy_matrix_part(H_copy, M * J - J + 1, K * M);
            break;
        case RU_code:
            cutted_H = copy_matrix(H_copy);
            break;
    }

    matrix G = create_G_from_H_matrix(cutted_H, columns_mdata);

    ldpc ldpc_object;
    ldpc_object.G = G;
    ldpc_object.H = H;
    ldpc_object.columns_mdata = columns_mdata;
    ldpc_object.n = G.columns;
    ldpc_object.k = G.rows;
    ldpc_object.r = H.rows;
    ldpc_object.systematic_r = H_copy.rows;
    ldpc_object.systematic_H = H_copy;
    ldpc_object.C = get_non_zero_column_data(H);
    ldpc_object.V = get_non_zero_column_data(transpose_matrix(H));

    free_matrix(cutted_H);

    return ldpc_object;
}

void free_ldpc(ldpc ldpc_object) {
    free_matrix(ldpc_object.G);
    free_matrix(ldpc_object.H);
	free(ldpc_object.columns_mdata.check_set);
	free(ldpc_object.columns_mdata.information_set);
}

matrix create_V_Gallager(int J, int K, int M) {
    int r = J * M;
    int n = K * M;
    int *x = (int *) malloc(n * sizeof(int));
    fill_with_permutation(x, n);
    matrix V = create_zero_matrix(r, K);
    int I = 0;
    int h, i, index;

    for (h = 0; h < J; h++) {
        int p = 0;
        for (i = 0; i < M; i++) {
            for (index = 0; index < K; index++) {
                V.body[I][index] = x[p + index];
            }
            p = p + K;
            I++;
        }
        fill_with_permutation(x, n);
    }


    free(x);
    return V;
}

matrix create_V_RU(int J, int K, int M) {
    int r = J * M;
    int n = K * M;
    int *A = (int *) malloc(n * J * sizeof(int));
    fill_with_permutation(A, n * J);
    int i, j;

    for (i = 0; i < n * J; i++) {
        A[i] = (A[i] - 1) / J;
    }

    matrix V = create_zero_matrix(J * M, K);

    for (j = 0; j < K; j++) {
        for (i = 0; i < r; i++) {
            V.body[i][j] = A[i + j * J * M];
        }
    }
    free(A);
    return V;
}


matrix create_H_rand(code_type type, int J, int K, int M) {
    int r = J * M;
    int n = K * M;
    matrix V;
    switch (type) {
        case Gallager:
            V = create_V_Gallager(J, K, M);
            break;
        case RU_code:
            V = create_V_RU(J, K, M);
            break;
    }

    int rw[r];
    int i, h;
    for (h = 0; h < r; h++) {
        rw[h] = 0;
        for (i = 0; i < K; i++) {
            rw[h] += V.body[h][i];
        }
    }

    matrix H = create_zero_matrix(r, n);

    for (i = 0; i < r; i++) {
        for (h = 0; h < K; h++) {
            H.body[i][V.body[i][h]] = 1;
        }
    }
    free_matrix(V);
    return H;

}

columns_metadata create_columns_metadata(int* check_set, int n, int k) {
	int R[n];
    int information_size = n;
    int i;
    for (i = 0; i < n; i++) {
        R[i] = 1;
    }
    for (i = 0; i < k; i++) {
        if (check_set[i] != -1) {
            R[check_set[i]] = 0;
            information_size--;
        }
    }

    int *information_set = (int *) malloc(information_size * sizeof(int));
    information_size = 0;
    for (i = 0; i < n; i++) {
        if (R[i] == 1) {
            information_set[information_size] = i;
            information_size++;
        }
    }

    check_set = (int*) realloc(check_set, (n - information_size) * sizeof(int));

    columns_metadata columns_mdata;
    columns_mdata.information_set = information_set;
    columns_mdata.check_set = check_set;
    columns_mdata.information_size = information_size;
    columns_mdata.check_size = n - information_size;

    return columns_mdata;
}

indices_of_nonzero_elements get_non_zero_column_data(matrix matrix_object) {
    indices_of_nonzero_elements result;
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
