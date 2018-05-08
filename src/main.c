#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "matrix.h"
#include "encoder.h"
#include "decoder.h"
#include "ldpc_generator.h"
#include "ldpc_tester.h"
#include  "subsidary_math.h"


#define TRUE !0
#define FALSE 0


void print_matrix(matrix M);
void print_ldpc(ldpc ldpc_object);


int main() {
    int J = 4, K = 8, M = 8;
    const int SUBMATRIX_SIZE = 67;
    //ldpc ldpc_object = create_ldpc(Gallager, J, K, M);
    int i = 0;
    int j = 0;
    int k = 0;
    //print_ldpc(ldpc_object);
    //SNR_interval SNR = {1., 5., 0.5};

    matrix weight_matrix = create_zero_matrix(2, 4);
    weight_matrix.body[0][0] = 3;
    weight_matrix.body[0][1] = 2;
    weight_matrix.body[0][2] = 3;
    weight_matrix.body[0][3] = 2;
    weight_matrix.body[1][0] = 0;
    weight_matrix.body[1][1] = 3;
    weight_matrix.body[1][2] = 2;
    weight_matrix.body[1][3] = 3;

    int num_of_weights = 0;
    int *num_of_weights_ptr = &num_of_weights;
    int ***polynomial_matrix = (int***)malloc(sizeof(int**) * weight_matrix.rows);
    for (i = 0; i < weight_matrix.rows; i++) {
        polynomial_matrix[i] = (int**)malloc(sizeof(int*) * weight_matrix.columns);
        for (j = 0; j < weight_matrix.columns; j++) {
            polynomial_matrix[i][j] = (int*)malloc(sizeof(int) * SUBMATRIX_SIZE);
            for (k = 0; k < SUBMATRIX_SIZE; k++) {
                polynomial_matrix[i][j][k] = 0;
            }
        }
    }

    weight_number_pair *w_n_pairs = get_weight_number_pairs(weight_matrix, num_of_weights_ptr);

    printf("%i\n", num_of_weights);
    printf("%i %i\n", w_n_pairs[0].weight, w_n_pairs[0].number);
    printf("%i %i\n", w_n_pairs[1].weight, w_n_pairs[1].number);
    get_polynomial_matrix(weight_matrix, SUBMATRIX_SIZE, num_of_weights, w_n_pairs, polynomial_matrix);
	printf("polynomial_matrix:\n");
    for (i = 0; i < weight_matrix.rows; i++) {
        for (j = 0; j < weight_matrix.columns; j++) {
            for (k = 0; k < SUBMATRIX_SIZE; k++) {
                printf("%i ", polynomial_matrix[i][j][k]);
            }
            printf("\n");
        }
    }
	get_polynomial_matrix_with_shift(polynomial_matrix, weight_matrix, SUBMATRIX_SIZE);
    printf("polynomial_matrix:\n");
    for (i = 0; i < weight_matrix.rows; i++) {
        for (j = 0; j < weight_matrix.columns; j++) {
            for (k = 0; k < SUBMATRIX_SIZE; k++) {
                printf("%i ", polynomial_matrix[i][j][k]);
            }
            printf("\n");
        }
    }

    print_matrix(create_H_matrix_use_polynomial_matrix_with_shifts(polynomial_matrix, weight_matrix, SUBMATRIX_SIZE));
    free(w_n_pairs);

    /*if (TRUE) {

    	FILE *file = fopen("decoding_simulation.txt", "w");
    	simulate_decoding(ldpc_object, SNR, file);
    	fclose(file);

	} else {

	    int k = ldpc_object.k, n = ldpc_object.n;
		int num_sigmas = (int) roundf((SNR.max - SNR.min) / SNR.step + 1);
		int i = 0;

		float *sigma_values = gen_sigma_values(SNR, (float) k / (float) n);
		printf("Sigma values:\n");
		for (i = 0; i < num_sigmas; i++) {
			printf("%.2f ", sigma_values[i]);
		}
	    printf("\n\n");

	    matrix U = create_random_matrix(1, k);
	    printf("U:\n");
	    print_matrix(U);
	    printf("\n");

	    matrix X = encode_message(ldpc_object, U);
	    printf("X:\n");
	    print_matrix(X);
	    printf("\n");

	    float *y = get_channel_output(X);
	    printf("Channel output:\n");
	    for (i = 0; i < n; i++) {
	    	printf("%.1f ", y[i]);
		}
		printf("\n\n");

		int SIGMA_INDEX = 2;
	    printf("Changes = %d\n", add_noise(y, n, sigma_values[SIGMA_INDEX]));
	    printf("Noised y:\n");
	    for (i = 0; i < n; i++) {
	    	printf("%.1f ", y[i]);
		}
		printf("\n\n");

	    normalize_message(y, n, sigma_values[SIGMA_INDEX] *
		                        sigma_values[SIGMA_INDEX]);
	    printf("Normalized y:\n");
	    for (i = 0; i < n; i++) {
	    	printf("%.1f ", y[i]);
		}
		printf("\n\n");

		for (i = 0; i < n; i++) {
	    	if (((y[i] < 0) ? 1 : 0) == X.body[0][i]) {
	    		printf(". ");
			} else {
				printf("# ");
			}
		}
		printf("\n\n");

		matrix *hard_solution = (matrix*)malloc(sizeof(matrix));
	    printf("Iterations = %d\n", flooding(ldpc_object, y, hard_solution));
	    printf("Hard solution:\n");
		print_matrix(*hard_solution);
	    printf("\n");

	    for (i = 0; i < n; i++) {
	    	if (hard_solution->body[0][i] == X.body[0][i]) {
	    		printf(". ");
			} else {
				printf("# ");
			}
		}
		printf("\n\n");

		printf("Syndrome:\n");
	    print_matrix(count_syndrome(ldpc_object, *hard_solution, TRUE));
	    printf("\n");
	}
*/
    system("pause");
    //free_ldpc(ldpc_object);

    return 0;
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


void print_ldpc(ldpc ldpc_object) {

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

	int i = 0;
    int j = 0;
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
        if (ldpc_object.C.element_length[i] == 0) printf("-");
        printf("\n");
    }
    printf("\n");

    printf("V =\n");
    for (i = 0; i < ldpc_object.H.rows; i++) {
        for (j = 0; j < ldpc_object.V.element_length[i]; j++) {
            printf("%d ", ldpc_object.V.element_data[i][j]);
        }
        if (ldpc_object.V.element_length[i] == 0) printf("-");
        printf("\n");
    }
}
