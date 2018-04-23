#include "matrix.h"
#include "ldpc_generator.h"

#ifndef DECODER_H_INCLUDED
#define DECODER_H_INCLUDED

float* map_sp(float y[], int length);
char check_syndrome(matrix hard, int r, indices_of_nonzero_elements V);
matrix get_hard_from_soft(float soft[], int length);
int flooding(ldpc ldpc_object, float *soft, matrix *hard_solution);
int decode_belief_propogandation(ldpc ldpc_object, float *y, matrix *hard_solution, char use_non_zero_data);
matrix count_syndrome(ldpc ldpc_object, matrix codedword, char use_non_zero_data);

#endif // DECODER_H_INCLUDED
