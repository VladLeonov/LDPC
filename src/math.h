#include "matrix.h"
#include "ldpc.h"

#ifndef MATH
#define MATH

int* gauss_elimination(matrix G);
int sum_row_elements(matrix G, int row_index);
void fill_with_permutation(int *x, int n);
matrix create_G_from_H_matrix(matrix H, columns_metadata columns_mdata);
int get_indexes_of_common_elements(int *arr_a, int *arr_b, int *result, int len_a, int len_b);
int decode_belief_propogandation(ldpc ldpc_object, float *y, matrix *hard_solution, char use_non_zero_data);
int flooding(ldpc ldpc_object, float *soft, matrix *hard_solution);

#endif
