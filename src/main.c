#include "ldpc.h"
#include "matrix.h"

int main() {

	/*matrix message = create_empty_matrix(1, 20);
	int i;
	for (i = 0; i < 20; i++) {
		message.body[0][i] = 'a' + i;
	}
	
	ldpc ldpc_object = create_ldpc(7, 4);
	matrix encoded = encode(ldpc_object, message);
	matrix decoded = decode(ldpc_object, encoded);
	
	for (i = 0; i < 20; i++) {
		printf("%c ", message.body[0][i]);
	}*/
	
	matrix M = create_zero_matrix(3, 5);
	M.body[0][0] = 1;
	M.body[0][1] = 1;
	M.body[0][3] = 1;
	M.body[1][1] = 1;
	M.body[1][2] = 1;
	M.body[1][4] = 1;
	M.body[2][0] = 1;
	M.body[2][2] = 1;
	M.body[2][3] = 1;
	
	
	printf("G = \n");
	print_matrix(M);
	print_matrix(copy_matrix(M));
	printf("\n");
	
	//create_systematic_view(M);
	
	system("pause");
	return 0;
}
