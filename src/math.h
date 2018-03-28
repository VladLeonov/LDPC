#ifndef MATH
#define MATH

int* gauss_elimination(matrix G);
int sum_rows(matrix G, int row_index);
void fill_with_permutation(int *x, int n);
matrix create_G_from_H_matrix(matrix H, columns_metadata columns_mdata);
int get_indexes_of_common_elements(int *arr_a, int *arr_b, int *result, int len_a, int len_b);
#endif
