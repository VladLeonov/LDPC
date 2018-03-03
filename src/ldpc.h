#include "matrix.h"

#ifndef LDPC
#define LDPC

typedef struct{
	matrix G, H;
	int n, k;
		
} ldpc;

matrix encode(ldpc ldpc_object, matrix message);
matrix decode(ldpc ldpc_object, matrix codeword);
matrix copy_matrix(matrix matrix_object);
ldpc create_ldpc(int n, int k);
void free_ldpc(ldpc ldpc_object);

#endif
