#include "encoder.h"


matrix encode_message(ldpc ldpc_object, matrix message) {
    int n = ldpc_object.n;
	int k = ldpc_object.k;
	int r = ldpc_object.systematic_r;
    matrix data = copy_matrix_part(message, 1, k);
    matrix codeword = create_zero_matrix(1, n);
    int* check_set = ldpc_object.columns_mdata.check_set;
    int* information_set = ldpc_object.columns_mdata.information_set;
    int i, j;

    for (i = 0; i < k; i++) {
    	codeword.body[0][information_set[i]] = data.body[0][i];
	}

	for (i = 0; i < r; i++) {
	    for (j = 0; j < k; j++) {
		    codeword.body[0][check_set[i]] ^= data.body[0][j] * ldpc_object.systematic_H.body[i][information_set[j]];
	    }
	}

	free_matrix(data);

	return codeword;
}
