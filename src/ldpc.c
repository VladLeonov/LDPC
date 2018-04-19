#include "ldpc.h"

#include "matrix.h"
#include "math.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define TRUE !0
#define FALSE 0

matrix encode(ldpc ldpc_object, matrix message, char use_non_zero_data) {
    int n = ldpc_object.n;
	int k = ldpc_object.k;
	int r = ldpc_object.r;
    matrix data = copy_matrix_part(message, 1, k);
    matrix codeword = create_zero_matrix(1, n);

    int* check_set = ldpc_object.columns_mdata.check_set;
    int* information_set = ldpc_object.columns_mdata.information_set;

    int i, j;
    for (i = 0; i < k; i++) {
    	codeword.body[0][information_set[i]] = data.body[0][i];
	}

	if (use_non_zero_data == FALSE) {

		for (i = 0; i < r; i++) {
		    for (j = 0; j < k; j++) {
			    codeword.body[0][check_set[i]] ^= data.body[0][j] * ldpc_object.systematic_H.body[i][information_set[j]];
		    }
	    }

	} else {

		int* data_indices = (int*)malloc(ldpc_object.columns_mdata.information_size * sizeof(int));
		int num_indices;
		for (i = 0; i < r; i++) {
			num_indices = get_indexes_of_common_elements(information_set, ldpc_object.V.element_data[i], data_indices,
			                                             ldpc_object.columns_mdata.information_size, ldpc_object.V.element_length[i]);

		    for (j = 0; j < num_indices; j++) {
			    codeword.body[0][check_set[i]] ^= data.body[0][data_indices[j]];
		    }
	    }

	}

	free_matrix(data);

	return codeword;
}

matrix count_syndrome(ldpc ldpc_object, matrix codedword, char use_non_zero_data) {

    matrix syndrome;
    if (use_non_zero_data == FALSE) {

    	matrix H_transposed = transpose_matrix(ldpc_object.H);
	    syndrome = multiply_matrices(codedword, H_transposed);
	    free_matrix(H_transposed);

    } else {

    	syndrome = create_zero_matrix(1, ldpc_object.r);
    	int i, j;
    	for (i = 0; i < ldpc_object.r; i++) {
    		for (j = 0; j < ldpc_object.V.element_length[i]; j++) {
    			syndrome.body[0][i] ^= codedword.body[0][ldpc_object.V.element_data[i][j]];
			}
		}

	}

    return syndrome;
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
    ldpc_object.r = cutted_H.rows;
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

void print_ldpc(ldpc ldpc_object) {
    int i = 0;
    int j = 0;

    printf("k = ");
    printf("%d\n\n", ldpc_object.k);

    printf("r = ");
    printf("%d\n\n", ldpc_object.r);

    printf("n = ");
    printf("%d\n\n", ldpc_object.n);

    printf("G =\n");
    print_matrix(ldpc_object.G);
    printf("\n");

    printf("H =\n");
    print_matrix(ldpc_object.H);
    printf("\n");

    columns_metadata columns_mdata = ldpc_object.columns_mdata;

    printf("Check set =\n");
    for (i = 0; i < columns_mdata.check_size; i++) {
        printf("%d ", columns_mdata.check_set[i]);
    }
    printf("\n\n");

    printf("Information set =\n");
    for (i = 0; i < columns_mdata.information_size; i++) {
        printf("%d ", columns_mdata.information_set[i]);
    }
    printf("\n\n");

    printf("C =\n");
    for (i = 0; i < ldpc_object.H.columns; i++) {
        for (j = 0; j < ldpc_object.C.element_length[i]; j++) {
            printf("%d ", ldpc_object.C.element_data[i][j]);
        }
        if (ldpc_object.C.element_length[i] == 0) {
        	printf("-");
		}
        printf("\n");
    }

    printf("\n");

    printf("V =\n");
    for (i = 0; i < ldpc_object.H.rows; i++) {
        for (j = 0; j < ldpc_object.V.element_length[i]; j++) {
            printf("%d ", ldpc_object.V.element_data[i][j]);
        }
        if (ldpc_object.V.element_length[i] == 0) {
        	printf("-");
		}
        printf("\n");
    }
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
        for (i = 0; i < J * M; i++) {
            V.body[i][j] = A[i + j * J * M];
        }
    }
    free(A);
    return V;
}

int get_max_element(int *array, int size) {
    int result = (int) pow(2, sizeof(int) * 8);
    int i;

    for (i = 0; i < size; i++) {
        if (result < array[i]) {
            result = array[i];
        }
    }

    return result;
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
    int i, j;
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

