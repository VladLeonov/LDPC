#include "ldpc_tester.h"

#include <stdlib.h>

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
