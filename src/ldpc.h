#include "matrix.h"

#ifndef LDPC
#define LDPC

typedef struct{
	matrix G, H;
	int n, k, check_size, information_size;
	int *check_set;
	int *information_set;
} ldpc;

matrix encode(ldpc ldpc_object, matrix message);
matrix count_syndrome(ldpc ldpc_object, matrix codedword);
ldpc create_ldpc(int n, int k);
void free_ldpc(ldpc ldpc_object);
void print_ldpc(ldpc ldpc_object);
matrix gen_LDPC_rand(int J, int K, int M);

#endif
