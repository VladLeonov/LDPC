#include "ldpc.h"
#include "matrix.h"

int main() {

	matrix message = create_empty_matrix(1, 20);
	int i;
	for (i = 0; i < 20; i++) {
		message.body[0][i] = 'a' + i;
	}
	
	ldpc ldpc_object = create_ldpc(7, 4);
	matrix encoded = encode(ldpc_object, message);
	matrix decoded = decode(ldpc_object, encoded);
	
	for (i = 0; i < 20; i++) {
		printf("%c ", message.body[0][i]);
	}
	system("pause");
	return 0;
}
