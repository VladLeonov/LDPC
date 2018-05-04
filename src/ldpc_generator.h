/**
    LDPC
    ldpc_generator.h
    Purpose: contains methods
    for generating LDPC code structure

    @author Danilyuk K.S.
    @version 25.04.2018
*/

#include "matrix.h"

#ifndef LDPC_GENERATOR_H_INCLUDED
#define LDPC_GENERATOR_H_INCLUDED


/**
    Enumeration contains types of LDPC codes
*/
typedef enum {
    Gallager,           //Gallager code
    RU_code             //RU_code
} code_type;

/**
    Struct contains data about information and check columns
*/
typedef struct {
    int check_size;             //size of array of check columns
    int information_size;       //size of array of information columns
    int *check_set;             //indices of check columns
    int *information_set;       //indices of information columns
} columns_metadata;

/**
    Struct contains data about indices of nonzero elements
*/
typedef struct {
    int **element_data;         //array which contains indices of nonzero elements
    int *element_length;        //array which contains number of non-zero
                                //elements in each column
} indices_of_nonzero_elements;

/**
    Struct contains data about LDPC code
*/
typedef struct {
    matrix G;                           //G matrix of LDPC code
    matrix H;                           //H matrix of LDPC code
    matrix systematic_H;                //Systematic H matrix of LDPC code
    int n;                              //Length of codeword
    int k;                              //Number of rows of G matrix
    int r;                              //Number of rows of H matrix
    int systematic_r;                   //Number of rows of systematic H matrix
    columns_metadata columns_mdata;     //data about information
                                        //and check columns of H matrix
    indices_of_nonzero_elements C;      //Array of numbers of nonzero
                                        //rows in each column
    indices_of_nonzero_elements V;      //Array of numbers of nonzero
                                        //columns in each row
} ldpc;

typedef struct {

    int weight;
    int number;

} weight_number_pair;

/**
    Creates LDPC code structure

    @param type The type of LDPC code
    @param J The number of ones per column
    @param K The number of ones per row
    @param M The size of submatrix
    @return The LDPC code structure
*/
ldpc create_ldpc(code_type type, int J, int K, int M);
/**
    Destroys LDPC code structure

    @param ldpc_object The structure of LDPC code
*/
void free_ldpc(ldpc ldpc_object);
/**
    Gets data about indices of nonero elements of H matrix

    @param matrix_object The matrix about which information is collected
    @return structure contaning data about nonero elements
*/
indices_of_nonzero_elements get_non_zero_column_data(matrix matrix_object);

 weight_number_pair* get_weight_number_pairs(matrix weight_matrix, int *num_of_weights_ptr);

#endif // LDPC_GENERATOR_H_INCLUDED
