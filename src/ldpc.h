#include "matrix.h"

#ifndef LDPC
#define LDPC

typedef struct{
	matrix G, H;
	int n, k;
		
} ldpc;

matrix encode(matrix G, matrix message);
matrix decode(matrix H, matrix codeword, int n, int k);
ldpc create_ldpc(int n, int k);
void free_ldpc(ldpc ldpc_object);

#endif
