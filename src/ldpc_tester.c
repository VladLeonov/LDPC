#include <stdlib.h>
#include <math.h>

#include "ldpc_tester.h"
#include "ldpc_generator.h"
#include "encoder.h"
#include "decoder.h"
#include "subsidary_math.h"


#define TRUE !0
#define FALSE 0


// Forms channel output from encoded message.
float* get_channel_output(matrix M) {
	float *channel_output = NULL;
	
	if (!is_void_matrix(M)) {
		channel_output = (float*) malloc(M.columns * sizeof(float));
		
		int i = 0;
		
		for (i = 0; i < M.columns; i++) {
			channel_output[i] = 2 * M.body[0][i] - 1;
		}
	}
	
	return channel_output;
}


//Calculates sigma value for each SRN value in interval.
void normalize_message(float *message, int length, float square_of_sigma) {
	int i = 0;
	
	for (i = 0; i < length; i++) {
		message[i] *= -2 / square_of_sigma;
	}
}


//Adds additive white Gaussian noise to message.
float* gen_sigma_values(SNR_interval SNRs, float R) {
	int array_size = (int) roundf((SNRs.max - SNRs.min) / SNRs.step + 1);
	float *sigma_values = (float*) malloc(array_size * sizeof(float));
	int i = 0;
	
	for (i = 0; i < array_size; i++) {
		sigma_values[i] = 
				sqrt(pow(10., -(SNRs.min + SNRs.step * i) / 10.) / 2. / R);
	}
	
	return sigma_values;
}


//Adds additive white Gaussian noise to message.
int add_noise(float *message, int length, float sigma) {
	int i = 0;
	float delta = 0;
	int changes = 0;
	
	for (i = 0; i < length; i++) {
		delta = sigma * randn();
		if ((message[i] + delta) * message[i] < 0) changes++;
		message[i] += delta;
	}
	
	return changes;
}


//Simulates the transmission of a message on the AWGN channel with encoding and decoding,
//write the simulation results (code characteristics) in a file.
void simulate_decoding(ldpc ldpc_object, SNR_interval SNRs,
                       FILE* output_file) {
    int k = ldpc_object.k, n = ldpc_object.n;
    float R = (float) k / (float) n;
    float *sigma_values = gen_sigma_values(SNRs, R);
    matrix *hard_solution = (matrix*)malloc(sizeof(matrix));
    int NEXP = 10000;
    int NERR = 100;
	float SNR = 0;
	int i = 0;
	
	fprintf(output_file, "SNR | FER for all messages | "
	                     "FER for incorrect messages | "
						 "average err in message | "
						 "average fixed err in message\n");
	
    for (SNR = SNRs.min, i = 0; SNR <= SNRs.max; SNR += SNRs.step, i++) {
    	printf("\nSNR = %f\n\n", SNR);
    	
        float FER = 0;
        float change_counter = 0;
        float changes_counter = 0;
        float fixed_errors = 0;
        int fixes = 0;
        int j = 0;
        
        for (j = 0; j < NEXP; j++) {
            matrix U = create_random_matrix(1, k);
            matrix X = encode_message(ldpc_object, U);
            float *y = get_channel_output(X);
            int changes = add_noise(y, n, sigma_values[i]);
			normalize_message(y, n, sigma_values[i] * sigma_values[i]);
            flooding(ldpc_object, y, hard_solution);
            
            if (changes > 0) {
            	change_counter += 1.0;
            	changes_counter += changes;
			}

            if ((j % 100 == 0) && (j > 0)) {
				printf("%02d%% iterations\n", j / 100);
			}

            if (compare_matrices(X, *hard_solution) == FALSE) {
            	FER += 1.0;
            	printf("%02d%% errors\n", (int) FER);
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

            if (FER >= NERR) {
            	j++;
                break;
            }
        }
        
        FER /= j;
        change_counter /= j;
        changes_counter /= j;

        fprintf(output_file, "%.1f %15.4f %24.4f %26.2f %26.2f\n",
		        SNR, FER, FER / change_counter,
				changes_counter / change_counter,
				fixed_errors / fixes);
    }
    
    free(hard_solution);
}
