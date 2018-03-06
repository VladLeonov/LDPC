#include "ldpc.h"
#include "matrix.h"

matrix encode(ldpc ldpc_object, matrix message) {
    return multiply_matrices(message, ldpc_object.G);
}

matrix count_syndrome(ldpc ldpc_object, matrix codedword) {
	
	matrix H = copy_matrix(ldpc_object.H);
    matrix H_transposed = transpose_matrix(H);
	matrix syndrome = multiply_matrices(codedword, H_transposed);
    
    return syndrome;
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
