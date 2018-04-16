#include "ldpc.h"
#include "matrix.h"
#include "math.h"
#include "ldpc_tester.h"

#include <stdlib.h>
#include <stdio.h>

#define TRUE !0
#define FALSE 0

int main() {

    int J = 2, K = 4, M = 3;//int J = 4, K = 8, M = 32;
    ldpc ldpc_object = create_ldpc(Gallager, J, K, M);
    
    SNR_interval SNR = {1., 5., 0.5};
    /*FILE *file = fopen("decoding_simulation.txt", "w");
    decoding_simulation(ldpc_object, SNR, file);
    fclose(file);*/
    
    matrix U, X;
    int k = ldpc_object.k, n = ldpc_object.n, i;
	float *y;
	matrix *hard_solution = (matrix*)malloc(sizeof(matrix));
	float *sigma_values = gen_sigma_values(SNR, (float) k / (float) n);
    
    printf("Information set:\n");
    for (i = 0; i < ldpc_object.columns_mdata.information_size; i++) {
    	printf("%d ", ldpc_object.columns_mdata.information_set[9]);
	}
	printf("\n\n");
    
    printf("U:\n");
    U = create_random_message(k);
    print_matrix(U);
    printf("\n");
    
    printf("X:\n");
    X = encode(ldpc_object, U, TRUE);
    print_matrix(X);
    printf("\n");
    
    printf("y:\n");
    y = normalize_vector(X, 1, -2);
    for (i = 0; i < n; i++) {
    	printf("%.1f ", y[i]);
	}
	printf("\n\n");
    
    printf("Noised y:\n");
    add_noise(y, n, sigma_values[0]);
    for (i = 0; i < n; i++) {
    	printf("%.1f ", y[i]);
	}
	printf("\n\n");
    
    
    printf("Iterations = %d\n", flooding(ldpc_object, y, hard_solution));
    printf("Hard solution:\n");
	print_matrix(*hard_solution);
    printf("\n");
	
    //int message_size = M * (K - J) + J - 1;
    
    system("pause");
    free_ldpc(ldpc_object);

    return 0;
}
