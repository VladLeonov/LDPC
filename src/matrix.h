#ifndef MATRIX
#define MATRIX

struct matrix {
	
	byte **body;
	int rows, columns;
	
	void free();
};

matrix create_zero_matrix(int rows, int columns);
matrix create_unit_matrix(int rows);
matrix create_sparse_matrix(int rows, int columns);

matrix multiply_matrix(matrix M1, matrix M2);
matrix transpose_matrix(matrix M);
matrix conbine_matrix(matrix M1, matrix M2);

#endif
