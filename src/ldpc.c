#include "ldpc.h"
#include "matrix.h"

matrix encode(ldpc ldpc_object, matrix message) {
    return multiply_matrices(message, ldpc_object.G);
}

matrix decode(ldpc ldpc_object, matrix codeword) {
    int decoded_rows = codeword.rows;
    int decoded_columns = ldpc_object.k;
    matrix decoded = create_empty_matrix(1, decoded_rows * decoded_columns);
    
    int i, j;
    for (i = 0; i < decoded_rows; i++) {
        for (j = 0; j < decoded_columns; j++) {
            decoded.body[i * decoded_columns + j] = codeword.body[i * decoded_columns + j];
        }
    }
    
    return decoded;
}

ldpc create_ldpc(int n, int k) {
    ldpc ldpc_object;
    ldpc_object.n = n;
    ldpc_object.k = k;
    
    matrix ik = create_unit_matrix(k);
    matrix ir = create_unit_matrix(n - k);
    matrix P = create_sparse_matrix(k, n - k);
    matrix Pt = transpose_matrix(P);
    ldpc_object.G = combine_matrices(ik, P);
    ldpc_object.H = combine_matrices(Pt, ir);
    
    free_matrix(ik);
    free_matrix(ir);
    free_matrix(P);
    free_matrix(Pt);
    
    return ldpc_object;
}

void free_ldpc(ldpc ldpc_object) {
    free_matrix(ldpc_object.G);
    free_matrix(ldpc_object.H);
	free(ldpc_object.check_set);
	free(ldpc_object.information_set);   
}

void print_ldpc(ldpc ldpc_object) {
    int i = 0;
    int j = 0;

    printf("k = ");
    printf("%d\n\n", ldpc_object.k);
    
    printf("n = ");
    printf("%d\n\n", ldpc_object.n);
    
    printf("G =\n");
    print_matrix(ldpc_object.G);
    printf("\n");
    
    printf("H =\n");
    print_matrix(ldpc_object.H);
    printf("\n");
    
    printf("Check set =\n");
    for (i = 0; i < ldpc_object.check_size; i++) {
        printf("%d ", ldpc_object.check_set[i]);
    }
    printf("\n\n");
    
    printf("Information set =\n");
    for (i = 0; i < ldpc_object.information_size; i++) {
        printf("%d ", ldpc_object.information_set[i]);
    }
    printf("\n\n");
}
