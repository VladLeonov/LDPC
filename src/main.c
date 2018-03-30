#include "ldpc.h"
#include "matrix.h"
#include "math.h"

#define TRUE !0
#define FALSE 0

int main() {

    int J = 2, K = 5, M = 4;
    ldpc ldpc_object = create_ldpc(J, K, M);
    //print_ldpc(ldpc_object);

    int message_size = M * (K - J) + J - 1;
	int message[1][message_size];
	int i;
	for (i = 0; i < message_size; i++) {
		message[0][i] = i % 2;
	}

	matrix encoded_message = encode(ldpc_object, array_to_matrix(1, message_size, message), TRUE);
    encoded_message.body[0][8] ^= 1;
    //matrix syndrome = count_syndrome(ldpc_object, encoded_message, FALSE);

    double *recieved_message = (double*) malloc(encoded_message.columns * sizeof(double));

	for (i = 0; i < encoded_message.columns; i++) {
        recieved_message[i] = encoded_message.body[0][i] == 0 ? -1 : 1;
	}

    matrix *hard_solution = (matrix*)malloc(sizeof(matrix));

    printf("Iterations = %d\n", decode_belief_propogandation(ldpc_object, recieved_message, hard_solution, FALSE));

	printf("encoded_message =\n");
    print_matrix(encoded_message);
    printf("\n");

    /*printf("syndrome =\n");
    print_matrix(syndrome);
    printf("\n");*/

    printf("hard_solution =\n");
    print_matrix(*hard_solution);
    printf("\n");

    system("pause");

    free_ldpc(ldpc_object);
    free_matrix(encoded_message);
    //free_matrix(syndrome);

    return 0;
}
