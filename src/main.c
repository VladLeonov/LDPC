#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "matrix.h"
#include "encoder.h"
#include "decoder.h"
#include "ldpc_generator.h"
#include "ldpc_tester.h"
#include  "subsidary_math.h"
#include "tanner_graph.h"


#define TRUE !0
#define FALSE 0


void print_matrix(matrix M);
void print_ldpc(ldpc ldpc_object);
void print_graph(graph G);


int main() {
    int J = 2, K = 4, M = 67;
    const int SUBMATRIX_SIZE = M;
    int i = 0;
    int j = 0;
    int k = 0;

	int array[2][4] = {{3,2,3,2}, {0,3,2,3}};
    matrix weight_matrix = array_to_matrix(J, K, array);

    /*int ***polynomial_matrix = (int***)malloc(sizeof(int**) * weight_matrix.rows);
    for (i = 0; i < weight_matrix.rows; i++) {
        polynomial_matrix[i] = (int**)malloc(sizeof(int*) * weight_matrix.columns);
        for (j = 0; j < weight_matrix.columns; j++) {
            polynomial_matrix[i][j] = (int*)malloc(sizeof(int) * SUBMATRIX_SIZE);
            for (k = 0; k < SUBMATRIX_SIZE; k++) {
                polynomial_matrix[i][j][k] = 0;
            }
        }
    }

	int num_of_weights = 0;
    weight_number_pair *w_n_pairs = get_weight_number_pairs(weight_matrix, &num_of_weights);
    get_polynomial_matrix(weight_matrix, SUBMATRIX_SIZE, num_of_weights, w_n_pairs, polynomial_matrix);
	get_polynomial_matrix_with_shift(polynomial_matrix, weight_matrix, SUBMATRIX_SIZE);
    matrix H = create_H_matrix_use_polynomial_matrix_with_shifts(polynomial_matrix, weight_matrix, SUBMATRIX_SIZE);
    free(w_n_pairs);
    */

    ldpc ldpc_object = create_ldpc(New_ldpc_code, 0, 0, SUBMATRIX_SIZE, weight_matrix);
    print_ldpc(ldpc_object);
    free_ldpc(ldpc_object);
    system("pause");

    return 0;
}


void print_matrix(matrix M) {
    int i, j;
    for (i = 0; i < M.rows; i++) {
        for (j = 0; j < M.columns; j++) {
            printf("%d ", M.body[i][j]);
        }
        printf("\n");
    }
}


void print_ldpc(ldpc ldpc_object) {

    printf("k = ");
    printf("%d\n\n", ldpc_object.k);

    printf("r = ");
    printf("%d\n\n", ldpc_object.r);

    printf("n = ");
    printf("%d\n\n", ldpc_object.n);

    printf("G =\n");
    print_matrix(ldpc_object.G);
    printf("\n");

    printf("H =\n");
    print_matrix(ldpc_object.H);
    printf("\n");

	int i = 0;
    int j = 0;
    columns_metadata columns_mdata = ldpc_object.columns_mdata;

    printf("Check set =\n");
    for (i = 0; i < columns_mdata.check_size; i++) {
        printf("%d ", columns_mdata.check_set[i]);
    }
    printf("\n\n");

    printf("Information set =\n");
    for (i = 0; i < columns_mdata.information_size; i++) {
        printf("%d ", columns_mdata.information_set[i]);
    }
    printf("\n\n");

    printf("C =\n");
    for (i = 0; i < ldpc_object.H.columns; i++) {
        for (j = 0; j < ldpc_object.C.element_length[i]; j++) {
            printf("%d ", ldpc_object.C.element_data[i][j]);
        }
        if (ldpc_object.C.element_length[i] == 0) printf("-");
        printf("\n");
    }
    printf("\n");

    printf("V =\n");
    for (i = 0; i < ldpc_object.H.rows; i++) {
        for (j = 0; j < ldpc_object.V.element_length[i]; j++) {
            printf("%d ", ldpc_object.V.element_data[i][j]);
        }
        if (ldpc_object.V.element_length[i] == 0) printf("-");
        printf("\n");
    }
}


void print_graph(graph G) {
	int i, j;
	for (i = 0; i < G.number_of_vertices; i++) {
		printf("%d(%d) - ", i, G.degrees_of_vertices[i]);
        for (j = 0; j < G.degrees_of_vertices[i]; j++) {
            printf("%d ", G.adjacency_list[i][j]);
        }
        printf("\n");
    }
}
