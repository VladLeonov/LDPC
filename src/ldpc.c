#include "ldpc.h"
#include "matrix.h"

matrix encode(matrix G, matrix message) {
	return multiply_matrices(message, G);
}

matrix decode(matrix H, matrix codeword, int n, int k) {
	int decoded_rows = codeword.rows;
	int decoded_columns = k;
	matrix decoded = create_empty_matrix(1, decoded_rows * decoded_columns);	
	
	int i, j;
	for (i = 0; i < decoded_rows; i++) {
	
		for (j = 0; j < decoded_columns; j++) {
			
			decoded.body[i * k + j] = codeword.body[i][j];
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
}

void free_ldpc(ldpc ldpc_object) {
	free_matrix(ldpc_object.G);
	free_matrix(ldpc_object.H);	
}
