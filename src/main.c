#include "ldpc.h"
#include "matrix.h"
#include "math.h"

#include <stdlib.h>
#include <stdio.h>

#define TRUE !0
#define FALSE 0

int main() {

    int J = 2, K = 5, M = 4;
    ldpc ldpc_object = create_ldpc(RU_code, J, K, M);

    int message_size = M * (K - J) + J - 1;
	int message[1][message_size];
	int i;
	for (i = 0; i < message_size; i++) {
		message[0][i] = i % 2;
	}

	matrix encoded_message = encode(ldpc_object, array_to_matrix(1, message_size, message), TRUE);

	printf("encoded_message =\n");
    print_matrix(encoded_message);
    printf("\n");

    float *recieved_message = (float*) malloc(encoded_message.columns * sizeof(float));

	for (i = 0; i < encoded_message.columns; i++) {
        recieved_message[i] = 2 - 4 * encoded_message.body[0][i];
	}

	recieved_message[8] *= -0.5;//-2

    matrix *hard_solution = (matrix*)malloc(sizeof(matrix));

    printf("Iterations = %d\n", decode_belief_propogandation(ldpc_object, recieved_message, hard_solution, FALSE));

    printf("hard_solution =\n");
    print_matrix(*hard_solution);

    for (i = 0; i < encoded_message.columns; i++) {
		printf(encoded_message.body[0][i] == hard_solution->body[0][i] ? ". " : "- ");
	}
	printf("\n\n");

	int *hard = (int*) malloc(encoded_message.columns * sizeof(int));
	printf("Iterations = %d\n", flooding(recieved_message, hard, ldpc_object));
	printf("hard =\n");
	for (i = 0; i < encoded_message.columns; i++) {
		printf("%d ", hard[i]);
	}
	printf("\n");
    for (i = 0; i < encoded_message.columns; i++) {
		printf(encoded_message.body[0][i] == hard[i] ? ". " : "- ");
	}

    system("pause");

    free_ldpc(ldpc_object);
    free_matrix(encoded_message);

    return 0;
}
