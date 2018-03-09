#ifndef MATH
#define MATH

int* gauss_elimination(matrix G);
ldpc create_systematic_view(matrix G_old, char is_H_entered);
int sum_rows(matrix G, int row_index);
void fill_with_permutation(int *x, int n);

#endif
