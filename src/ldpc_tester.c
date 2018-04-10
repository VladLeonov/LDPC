#include "ldpc_tester.h"

#include <stdlib.h>
#include <math.h>

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

