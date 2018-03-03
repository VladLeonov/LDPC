#include "matrix.h"

#ifndef LDPC
#define LDPC

typedef struct{
	matrix G, H;
	int n, k, check_size;
	int *check_set;
	int *information_set;
		
} ldpc;

matrix encode(ldpc ldpc_object, matrix message);
matrix decode(ldpc ldpc_object, matrix codeword);
matrix copy_matrix(matrix matrix_object);
ldpc create_ldpc(int n, int k);
void free_ldpc(ldpc ldpc_object);
void print_ldpc(ldpc ldpc_object);

#endif
