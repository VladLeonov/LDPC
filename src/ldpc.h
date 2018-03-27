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
	int n, k, r;
	columns_metadata columns_mdata;
	non_zero_data C, V;
} ldpc;

matrix encode(ldpc ldpc_object, matrix message, char use_non_zero_data);
matrix count_syndrome(ldpc ldpc_object, matrix codedword);
ldpc create_ldpc(int J, int K, int M);
void free_ldpc(ldpc ldpc_object);
void print_ldpc(ldpc ldpc_object);
matrix create_H_rand(int J, int K, int M);
columns_metadata create_columns_metadata(int* information_set, int n, int k);

#endif
