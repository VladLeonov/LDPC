#include "matrix.h"
#include "ldpc_generator.h"


#ifndef DECODER_H_INCLUDED
#define DECODER_H_INCLUDED

int flooding(ldpc ldpc_object, float *soft, matrix *hard_solution);
matrix count_syndrome(ldpc ldpc_object, matrix codeword, char use_non_zero_data);

#endif // DECODER_H_INCLUDED
