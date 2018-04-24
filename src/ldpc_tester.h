#include "stdio.h"
#include "matrix.h"
#include "ldpc_generator.h"

#ifndef LDPC_TESTER_CORRECT_VERSION_H_INCLUDED
#define LDPC_TESTER_CORRECT_VERSION_H_INCLUDED

typedef struct{
	float min, max, step;
} SNR_interval;

float* get_channel_output(matrix M);
void normalize_message(float *message, int length, float square_of_sigma);
float* gen_sigma_values(SNR_interval SNR, float R);
int add_noise(float *message, int length, float sigma);
void simulate_decoding(ldpc ldpc_object, SNR_interval SNRs, FILE* output_file);

#endif // LDPC_TESTER_CORRECT_VERSION_H_INCLUDED
