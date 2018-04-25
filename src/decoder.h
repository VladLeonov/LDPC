#include "matrix.h"
#include "ldpc_generator.h"


#ifndef DECODER_H_INCLUDED
#define DECODER_H_INCLUDED

/**
    LDPC
    encoder.h
    Purpose: contains methods
    for decoding message by
    using LDPC code structure

    @author Danilyuk K.S.
    @version 25.04.2018
*/

/**
    Decodes accepted sequence by belief propogandation method

    @param ldpc_object The structure of LDPC code
    @param soft The pointer to array of float numbers
    which contains soft solution
    @param hard_solution The pointer to matrix which contains hard solution
    @return The number of iterations performed
*/
int flooding(ldpc ldpc_object, float *soft, matrix *hard_solution);
/**
    Counts syndrome for codeword which is result
    of coding message by LDPC code

    @param The ldpc_object structure of LDPC code
    @param The codeword pointer to matrix which contains
    sequence of ones and zeros
    @param The use_non_zero_data flag of using data about non-zero elements
    of H matrix of LDPC code
*/
matrix count_syndrome(ldpc ldpc_object, matrix codeword, char use_non_zero_data);

#endif // DECODER_H_INCLUDED
