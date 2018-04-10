#include "matrix.h"

#ifndef LDPC
#define LDPC

typedef enum code_type {
    Gallager = 0,
    RU_code = 1
} code_type;

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
matrix count_syndrome(ldpc ldpc_object, matrix codedword, char use_non_zero_data);
ldpc create_ldpc(code_type type, int J, int K, int M);
void free_ldpc(ldpc ldpc_object);
void print_ldpc(ldpc ldpc_object);
matrix create_H_rand(code_type type, int J, int K, int M);
columns_metadata create_columns_metadata(int* information_set, int n, int k);
matrix create_random_message(int length);

#endif
