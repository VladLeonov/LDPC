/**
    LDPC
    ldpc_tester.h
    Purpose: Simulates the transmission of a message on the AWGN channel
	in order to obtain the dependence of the error probability when decoding
	the message from the signal-to-noise ratio.

    @author Leonov V.R.
    @version 24.04.18
*/


#include "stdio.h"
#include "matrix.h"
#include "ldpc_generator.h"


#ifndef LDPC_TESTER_CORRECT_VERSION_H_INCLUDED
#define LDPC_TESTER_CORRECT_VERSION_H_INCLUDED


//Ñompactly stores a set of equally spaced values of signal-to-noise ratio.
typedef struct{
	float min;		//Min value.
	float max;		//Max value.
	float step;		//Space between values.
} SNR_interval;

/**
    Forms channel output from encoded message.

    @param M Encoded message.
    @return Channel output.
*/
float* get_channel_output(matrix M);

/**
    Adduce message to the normalized form.

    @param message Encoded message.
    @param length Message length.
    @param square_of_sigma Ñoefficient of normalization.
*/
void normalize_message(float *message, int length, float square_of_sigma);

/**
    Calculates sigma value for each SRN value in interval.

    @param SNRs Interval of SRN for which it is necessary 
	to calculate the sigma.
    @param R Code speed.
    @return Array of sigma values.
*/
float* gen_sigma_values(SNR_interval SNRs, float R);

/**
    Adds additive white Gaussian noise to message.

    @param message Encoded message.
    @param length Message length.
    @param sigma Ñoefficient of noise.
*/
int add_noise(float *message, int length, float sigma);

/**
    Simulates the transmission of a message on the AWGN channel with encoding and decoding,
	write the simulation results (code characteristics) in a file.

    @param ldpc_object Ñode used for encoding and decoding.
    @param SNRs Interval of SRN for which it is necessary 
	to calculate the code characteristics.
    @param output_file File for writing results.
*/
void simulate_decoding(ldpc ldpc_object, SNR_interval SNRs, 
                       FILE* output_file);

#endif // LDPC_TESTER_CORRECT_VERSION_H_INCLUDED
