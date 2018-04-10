#include "matrix.h"

#ifndef LDPC_TESTER
#define LDPC_TESTER

typedef struct{
	float min, max, step;
} SNR_interval;

matrix create_random_message(int length);
float* normalize_vector(matrix M, float shift, float factor);
float* gen_sigma_values(SNR_interval SNR, float R);
float randn();
void add_noise(float *message, int length, float sigma);

#endif
