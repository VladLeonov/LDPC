#include "ldpc_tester.h"

#include <stdlib.h>
#include <math.h>

#include "math.h"

#define TRUE !0
#define FALSE 0

matrix create_random_message(int length) {
    int i;
    matrix message = create_zero_matrix(1, length);
    for (i = 0; i < length; i++) {
        message.body[0][i] = rand() % 2;
    }
    return message;
}

float* normalize_vector(matrix M, float shift, float factor) {
	float *result = NULL;
	if (!is_void_matrix(M)) {
		result = (float*) malloc(M.columns * sizeof(float));
		int i;
		for (i = 0; i < M.columns; i++) {
			result[i] = M.body[0][i] * factor + shift;
		}
	}
	return result;
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

float randn() {
	float U1, U2, W, mult;
	static float X1, X2;
	static int call = 0;

	if (call == 1) {
		call = !call;
		return X2;
	}

	do {
		U1 = -1 + ((float) rand () / RAND_MAX) * 2;
		U2 = -1 + ((float) rand () / RAND_MAX) * 2;
		W = pow (U1, 2) + pow (U2, 2);
	} while ((W >= 1) || (W == 0));

	mult = sqrt ((-2 * log (W)) / W);
	X1 = U1 * mult;
	X2 = U2 * mult;

	call = !call;

	return X1;
}

void add_noise(float *message, int length, float sigma) {
	int i;
	for (i = 0; i < length; i++) {
		message[i] += sigma * randn();
	}
}

void decoding_simulation(ldpc ldpc_object, SNR_interval SNRs, FILE* output_file) {
	int NEXP = 10000;
    int NERR = 100;
	float PER;

    int k = ldpc_object.k, n = ldpc_object.n;
    float R = (float) k / (float) n;
    float *sigma_values = gen_sigma_values(SNRs, R);

	float SNR;
	matrix U, X;
	float *y;
	int i, j;
	matrix *hard_solution = (matrix*)malloc(sizeof(matrix));
    for (SNR = SNRs.min, i = 0; SNR <= SNRs.max; SNR += SNRs.step, i++) {
        PER = 0;
        for (j = 0; j < NEXP; j++) {
            U = create_random_message(k);
            X = encode(ldpc_object, U, TRUE);
            y = normalize_vector(X, 1, -2);
            add_noise(y, n, sigma_values[i]);
            flooding(ldpc_object, y, hard_solution);

            if (compare_matrices(X, *hard_solution) == FALSE) {
            	PER += 1.0;
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
        fprintf(output_file, "%f %f\n", SNR, PER);
    }
    free(hard_solution);
}
