#include "matrix.h"

#ifndef LDPC_GENERATOR_H_INCLUDED
#define LDPC_GENERATOR_H_INCLUDED

typedef enum {
    Gallager,
    RU_code
} code_type;

typedef struct {
    int check_size, information_size;
    int *check_set;
    int *information_set;
} columns_metadata;

typedef struct{
    int **element_data;
    int *element_length;
} indices_of_nonzero_elements;

typedef struct {
    matrix G, H, systematic_H;
    int n, k, r, systematic_r;
    columns_metadata columns_mdata;
    indices_of_nonzero_elements C, V;
} ldpc;

int* gauss_elimination();
matrix create_G_from_H_matrix();
ldpc create_ldpc();
void free_ldpc();
matrix create_V_Gallager();
matrix create_V_RU();
matrix create_H_rand();
columns_metadata create_columns_metadata();
indices_of_nonzero_elements get_non_zero_column_data();

#endif // LDPC_GENERATOR_H_INCLUDED
