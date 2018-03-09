#include "matrix.h"

#ifndef LDPC
#define LDPC

typedef struct{
	int check_size, information_size;
	int *check_set;
	int *information_set;
} columns_metadata;

typedef struct{
	matrix G, H;
	int n, k;
	columns_metadata columns_mdata;
} ldpc;

matrix encode(ldpc ldpc_object, matrix message);
matrix count_syndrome(ldpc ldpc_object, matrix codedword);
ldpc create_ldpc(int n, int k);
void free_ldpc(ldpc ldpc_object);
void print_ldpc(ldpc ldpc_object);
matrix create_H_rand(int J, int K, int M);

#endif
