#include "matrix.h"
#include "ldpc_generator.h"


#ifndef ENCODER_H_INCLUDED
#define ENCODER_H_INCLUDED

/**
    LDPC
    encoder.h
    Purpose: contains methods
    for encoding message by
    using LDPC code structure

    @author Danilyuk K.S.
    @version 25.04.2018
*/

/**
    Encodes message by LDPC code

    @param ldpc_object The structure of LDPC code
    @param message The matrix containing single row,
    which ones and zeros sequence
    @return The matrix containing single row, which
    is codeword of LDPC code
*/
matrix encode_message(ldpc ldpc_object, matrix message);

#endif // ENCODER_H_INCLUDED
