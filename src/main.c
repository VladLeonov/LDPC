#include "ldpc.h"
#include "matrix.h"
#include "math.h"
#include "ldpc_tester.h"

#include <stdlib.h>
#include <stdio.h>

#define TRUE !0
#define FALSE 0

int main() {

    int J = 4, K = 8, M = 32;
    ldpc ldpc_object = create_ldpc(Gallager, J, K, M);
    //print_ldpc(ldpc_object);

    SNR_interval SNR = {1., 9., 1.};

    if (TRUE) {

    	FILE *file = fopen("decoding_simulation.txt", "w");
    	decoding_simulation(ldpc_object, SNR, file);
    	fclose(file);

	} else {

	    matrix U, X;
	    int k = ldpc_object.k, n = ldpc_object.n, i;
	    int sigma_index = 2;
		float *y;
		matrix *hard_solution = (matrix*)malloc(sizeof(matrix));
		float *sigma_values = gen_sigma_values(SNR, (float) k / (float) n);
		int num_sigmas = (int) roundf((SNR.max - SNR.min) / SNR.step + 1);

		printf("Sigma values:\n");
		for (i = 0; i < num_sigmas; i++) {
			printf("%.2f ", sigma_values[i]);
		}
	    printf("\n\n");

	    U = create_random_message(k);
	    printf("U:\n");
	    print_matrix(U);
	    printf("\n");

	    X = encode(ldpc_object, U);
	    printf("X:\n");
	    print_matrix(X);
	    printf("\n");

	    y = get_channel_output(X);
	    printf("Channel output:\n");
	    for (i = 0; i < n; i++) {
	    	printf("%.1f ", y[i]);
		}
		printf("\n\n");

	    printf("Changes = %d\n", add_noise(y, n, sigma_values[sigma_index]));
	    printf("Noised y:\n");
	    for (i = 0; i < n; i++) {
	    	printf("%.1f ", y[i]);
		}
		printf("\n\n");

	    normalize_vector(y, n, sigma_values[sigma_index] * sigma_values[sigma_index]);
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

	    printf("Iterations = %d\n", flooding(ldpc_object, y, hard_solution));
	    //printf("Iterations = %d\n", decode_belief_propogandation(ldpc_object, y, hard_solution, TRUE));
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

    //int message_size = M * (K - J) + J - 1;

    system("pause");
    free_ldpc(ldpc_object);

    return 0;
}
