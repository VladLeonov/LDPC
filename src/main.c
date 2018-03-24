#include "ldpc.h"
#include "matrix.h"
#include "math.h"

int main() {

    int J = 2, K = 5, M = 4;
    ldpc ldpc_object = create_ldpc(J, K, M);
    print_ldpc(ldpc_object);
    printf("\n");
    
    int message_size = M * (K - J) + J - 1;
	int message[1][message_size];
	int i;
	for (i = 0; i < message_size; i++) {
		message[0][i] = i % 2;
	}
	
	matrix encoded_message = encode(ldpc_object, array_to_matrix(1, message_size, message));
	printf("Encoded message =\n");
	print_matrix(encoded_message);
	printf("\n");
	
	matrix syndrome = count_syndrome(ldpc_object, encoded_message);
	printf("Syndrome =\n");
	print_matrix(syndrome);
	printf("\n");
	
    system("pause");
    return 0;
}
