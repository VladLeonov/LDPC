#include "ldpc.h"
#include "matrix.h"

matrix encode(ldpc ldpc_object, matrix message) {
    return multiply_matrices(message, ldpc_object.G);
}

matrix count_syndrome(ldpc ldpc_object, matrix codedword) {
	
    matrix H_transposed = transpose_matrix(ldpc_object.H);
	matrix syndrome = multiply_matrices(codedword, H_transposed);
	free_matrix(H_transposed);
    
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
	free(ldpc_object.columns_mdata.check_set);
	free(ldpc_object.columns_mdata.information_set);   
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
}

matrix create_H_rand(int J, int K, int M) {
	int r = J * M;
	int n = K * M;
	int *x = (int *) malloc(n * sizeof(int));
	fill_with_permutation(x, n);
	matrix V = create_zero_matrix(r, K);
	int I = 0;
	int h, i, index;
	
	for (h = 0; h < J; h++) {
		int p = 0;
		for (i = 0; i < M; i++) {
			for (index = 0; index < K; index++) {
				V.body[I][index] = x[p + index];
			}
			p = p + K;
			I++;
		}
		fill_with_permutation(x, n);
	}
	
	int rw[r];
	
	for (h = 0; h < r; h++) {
		rw[h] = 0;
		for (i = 0; i < K; i++) {
			rw[h] += V.body[h][i];
		}
	}
	
	matrix H = create_zero_matrix(r, n);
	
	for (i = 0; i < r; i++) {
		for (h = 0; h < K; h++) {
			H.body[i][V.body[i][h]] = 1;
		}
	}
	
	free_matrix(V);
	
	return H;
}

columns_metadata create_columns_metadata(int* information_set, int n, int k) {
	int R[n];
    int check_size = n;
    int i, j;
    for (i = 0; i < n; i++) {
        R[i] = 1;
    }
    for (i = 0; i < k; i++) {
        if (information_set[i] != -1) {
            R[information_set[i]] = 0;
            check_size--;
        }
    }
    
    int *check_set = (int *) malloc(check_size);
    check_size = 0;
    for (i = 0; i < n; i++) {
        if (R[i] == 1) {
            check_set[check_size] = i;
            check_size++;
        }
    }
    
    columns_metadata columns_mdata;
    columns_mdata.information_set = information_set;
    columns_mdata.check_set = check_set;
    columns_mdata.information_size = n - check_size;
    columns_mdata.check_size = check_size;
    
    return columns_mdata;
}
