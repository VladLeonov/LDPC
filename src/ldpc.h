#include "matrix.h"

#ifndef LDPC
#define LDPC

struct ldpc {
	
	matrix G, H;
	int n, k;
	
	matrix encode(matrix message);
	matric decode(matrix codeword);
	void free();
};

ldpc create_ldpc(int n, int k);

#endif
