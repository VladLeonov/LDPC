#include <stdlib.h>
#include <math.h>
#include "ldpc_tester.h"
#include "ldpc_generator.h"
#include "encoder.h"
#include "decoder.h"
#include "subsidary_math.h"

#define TRUE !0
#define FALSE 0

float* get_channel_output(matrix M) {
	float *result = NULL;
	if (!is_void_matrix(M)) {
		result = (float*) malloc(M.columns * sizeof(float));
		int i;
		for (i = 0; i < M.columns; i++) {
			result[i] = 2 * M.body[0][i] - 1;
		}
	}
	return result;
}

void normalize_message(float *message, int length, float square_of_sigma) {
	int i;
	for (i = 0; i < length; i++) {
		message[i] *= -2 / square_of_sigma;
	}
}

float* gen_sigma_values(SNR_interval SNR, float R) {
	int array_size = (int) roundf((SNR.max - SNR.min) / SNR.step + 1);
	float *result = (float*) malloc(array_size * sizeof(float));
	int i;
	for (i = 0; i < array_size; i++) {
		result[i] = sqrt(pow(10., -(SNR.min + SNR.step * i) / 10.) / 2. / R);
	}
	return result;
}

int add_noise(float *message, int length, float sigma) {
	int i;
	float delta;
	int changes = 0;
	for (i = 0; i < length; i++) {
		delta = sigma * randn();
		if ((message[i] + delta) * message[i] < 0) {
			changes++;
		}
		message[i] += delta;
	}
	return changes;
}

void decoding_simulation(ldpc ldpc_object, SNR_interval SNRs, FILE* output_file) {

	int NEXP = 10000;
    int NERR = 100;
	float PER, change_counter, changes_counter, fixed_errors;

    int k = ldpc_object.k, n = ldpc_object.n;
    float R = (float) k / (float) n;
    float *sigma_values = gen_sigma_values(SNRs, R);

	float SNR;
	matrix U, X;
	float *y;
	int i, j;
	int changes, fixes;
	matrix *hard_solution = (matrix*)malloc(sizeof(matrix));
	fprintf(output_file, "SNR | PER for all messages | PER for incorrect messages | average err in message | average fixed err in message\n");
    for (SNR = SNRs.min, i = 0; SNR <= SNRs.max; SNR += SNRs.step, i++) {
    	printf("\nSNR = %f\n\n", SNR);
        PER = 0;
        change_counter = 0;
        changes_counter = 0;
        fixed_errors = 0;
        fixes = 0;
        for (j = 0; j < NEXP; j++) {
            U = create_random_matrix(1, k);
            X = encode(ldpc_object, U);
            y = get_channel_output(X);
            changes = add_noise(y, n, sigma_values[i]);
            if (changes > 0) {
            	change_counter += 1.0;
            	changes_counter += changes;
			}
			normalize_message(y, n, sigma_values[i] * sigma_values[i]);
            flooding(ldpc_object, y, hard_solution);
            //decode_belief_propogandation(ldpc_object, y, hard_solution, TRUE);

            if ((j % 100 == 0) && (j > 0)) {
				printf("%02d%% iterations\n", j / 100);
			}

            if (compare_matrices(X, *hard_solution) == FALSE) {
            	PER += 1.0;
            	printf("%02d%% errors\n", (int) PER);
			} else {
				if (changes > 0) {
					fixed_errors += changes;
					fixes++;
				}
			}

			free_matrix(U);
			free_matrix(X);
			free(y);
			free_matrix(*hard_solution);

            if (PER >= NERR) {
            	j++;
                break;
            }
        }
        PER /= j;
        change_counter /= j;
        changes_counter /= j;
        printf("%.1f %15.4f %24.4f %26.2f %26.2f\n", SNR, PER, PER / change_counter, changes_counter / change_counter, fixed_errors / fixes);


        fprintf(output_file, "%.1f %15.4f %24.4f %26.2f %26.2f\n", SNR, PER, PER / change_counter, changes_counter / change_counter, fixed_errors / fixes);
    }
    free(hard_solution);
}
