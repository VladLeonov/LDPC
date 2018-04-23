#include "matrix.h"

#include "ldpc.h"

#include <stdio.h>

#ifndef LDPC_TESTER
#define LDPC_TESTER

typedef struct{
	float min, max, step;
} SNR_interval;

matrix create_random_matrix(int rows, int columns);
float* get_channel_output(matrix M);
void normalize_message(float *message, int length, float square_of_sigma);
float* gen_sigma_values(SNR_interval SNR, float R);
float randn();
int add_noise(float *message, int length, float sigma);
void decoding_simulation(ldpc ldpc_object, SNR_interval SNRs, FILE* output_file);

#endif
