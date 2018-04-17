#include "ldpc.h"
#include "matrix.h"
#include "math.h"
#include "ldpc_tester.h"

#include <stdlib.h>
#include <stdio.h>

#define TRUE !0
#define FALSE 0

int main() {

    int J = 4, K = 8, M = 8;//J = 2, K = 4, M = 3;
    ldpc ldpc_object = create_ldpc(RU_code, J, K, M);
    
    SNR_interval SNR = {1., 9., 1.};
    /*FILE *file = fopen("decoding_simulation.txt", "w");
    decoding_simulation(ldpc_object, SNR, file);
    fclose(file);*/
    
    matrix U, X;
    int k = ldpc_object.k, n = ldpc_object.n, i;
    int sigma_index = 6;
	float *y;
	matrix *hard_solution = (matrix*)malloc(sizeof(matrix));
	float *sigma_values = gen_sigma_values(SNR, (float) k / (float) n);
    
    U = create_random_message(k);
    printf("U:\n");
    print_matrix(U);
    printf("\n");
    
    X = encode(ldpc_object, U, TRUE);
    printf("X:\n");
    print_matrix(X);
    printf("\n");
    
    y = get_channel_output(X);
    printf("Channel output:\n");
    for (i = 0; i < n; i++) {
    	printf("%.1f ", y[i]);
	}
	printf("\n\n");
    
    printf("Changes = %d\n", add_noise(y, n, sigma_values[sigma_index]));
    printf("Noised y:\n");
    for (i = 0; i < n; i++) {
    	printf("%.1f ", y[i]);
	}
	printf("\n\n");
	
    normalize_vector(y, n, sigma_values[sigma_index] * sigma_values[sigma_index]);
    printf("normalized y:\n");
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
