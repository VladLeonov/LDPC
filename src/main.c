#include "ldpc.h"
#include "matrix.h"
#include "math.h"

int main() {

    int J = 2, K = 5, M = 4;
    matrix H = create_H_rand(J, K, M);
	ldpc ldpc_object = create_systematic_view(H, 1);
	
	int message_size = M * (K - J);
	int message[1][message_size];
	int i;
	for (i = 0; i < message_size; i++) {
		message[0][i] = i % 2;
	}
	
	matrix encoded_message = encode(ldpc_object, array_to_matrix(1, message_size, message));
	matrix syndrome = count_syndrome(ldpc_object, encoded_message);
	
	printf("J = %d, K = %d, M = %d\n\n", J, K, M);
	printf("Random H =\n");
	print_matrix(H);
	printf("\n");
	print_ldpc(ldpc_object);
	printf("Syndrome =\n");
	print_matrix(syndrome);
	
    system("pause");
    
    free_matrix(H);
    free_matrix(encoded_message);
    free_matrix(syndrome);
    return 0;
}
