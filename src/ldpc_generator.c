#include <stdlib.h>

#include "ldpc_generator.h"
#include "matrix.h"
#include "subsidary_math.h"

#define TRUE !0
#define FALSE 0

/**
    Transform linearly dependent rows to zero rows
*/
int* perform_gauss_elimination(matrix G) {
	int k = G.rows;
	int n = G.columns;
    int* information_set = (int*) malloc(k * sizeof(int));
    int i, j;
    int column_number = 0;

    for (i = 0; i < k; i++){
        information_set[i] = -1;
    }

    i = 0;

    while (i < k) {

        //throw away all-zero row
        while (calculate_sum_row_elements(G, i) == 0) {
            for (j = i; j < (k - 1); j++) {
                G.body[j] = G.body[j + 1];
            }
            k--;
            information_set[k] = -1;
            if (i >= k) break;
        }

        if (i >= k) break;

        //find 1 in the column
        int b = 0;
        j = i - 1;

        while ((b == 0) && (j < k - 1)) {
            j++;
            b = G.body[j][column_number];
        }

        if (b == 1) {
            information_set[i] = column_number;

            // remove other ones from this column
            int o, u;
            for (o = 0; o < k; o++) {
                if (o == j) continue;
                if (G.body[o][column_number]==1) {
                    for (u = 0; u < n; u++) {
                        G.body[o][u] ^= G.body[j][u];
                    }
                }
            }

            // replace rows i and j
            int *buff = G.body[i];
            G.body[i] = G.body[j];
            G.body[j] = buff;

            i++;
        }

        column_number++;  // keep the same row, another column
    }

    for (i = k; i < G.rows; i++) {
    	G.body[i] = calloc(G.columns, sizeof(int));
	}

    return information_set;
}

/**
    Creates G matrix from H matrix
*/
matrix create_G_from_H_matrix(matrix H, columns_metadata columns_mdata) {
    matrix P = create_empty_matrix(H.rows, columns_mdata.information_size);
    int i = 0, j = 0;

    for (j = 0; j < columns_mdata.information_size; j++) {
        for (i = 0; i < H.rows; i++) {
            P.body[i][j] = H.body[i][columns_mdata.information_set[j]];
        }
    }

    matrix PT = transpose_matrix(P);
    matrix G = create_zero_matrix(H.columns - H.rows, H.columns);
    matrix U = create_unit_matrix(G.rows);

    for (j = 0; j < columns_mdata.check_size; j++) {
        for (i = 0; i < G.rows; i++) {
            G.body[i][columns_mdata.check_set[j]] = PT.body[i][j];
        }
    }

    for (j = 0; j < columns_mdata.information_size; j++) {
        for (i = 0; i < G.rows; i++) {
            G.body[i][columns_mdata.information_set[j]] = U.body[i][j];
        }
    }

    free_matrix(P);
    free_matrix(PT);
    free_matrix(U);

    return G;
}

/**
    Prepares data for creation H matrix of Gallager's code
*/
matrix create_V_Gallager(int J, int K, int M) {
    int r = J * M;
    int n = K * M;
    int *x = (int *) malloc(n * sizeof(int));
    matrix V = create_zero_matrix(r, K);
    int I = 0;
    int h = 0, i = 0, index = 0;

    fill_with_permutation(x, n);

    for (h = 0; h < J; h++) {
        int p = 0;
        for (i = 0; i < M; i++) {
            for (index = 0; index < K; index++) {
                V.body[I][index] = x[p + index];
            }
            p = p + K;
            I++;
        }
        fill_with_permutation(x, n);
    }

    free(x);

    return V;
}

/**
    Prepares data for creation H matrix of RU code
*/
matrix create_V_RU(int J, int K, int M) {
    int r = J * M;
    int n = K * M;
    int *A = (int *) malloc(n * J * sizeof(int));
    int i = 0, j = 0;
    matrix V = create_zero_matrix(r, K);

    fill_with_permutation(A, n * J);

    for (i = 0; i < n * J; i++) {
        A[i] = (A[i] - 1) / J;
    }

    for (j = 0; j < K; j++) {
        for (i = 0; i < r; i++) {
            V.body[i][j] = A[i + j * J * M];
        }
    }

    free(A);

    return V;
}

/**
    Creates random H matrix
*/
matrix create_H_rand(code_type type, int J, int K, int M) {

    int r = J * M;
    int n = K * M;
    matrix V;

    switch (type)
    {
        case GALLAGER:
            V = create_V_Gallager(J, K, M);
            break;
        case RU_CODE:
            V = create_V_RU(J, K, M);
            break;
        default:
            V = create_random_matrix(r, K);
    }

    int rw[r];
    int i, h;

    for (h = 0; h < r; h++) {
        rw[h] = 0;
        for (i = 0; i < K; i++) {
            rw[h] += V.body[h][i];
        }
    }

    matrix H = create_zero_matrix(r, n);

    for (i = 0; i < r; i++) {
        for (h = 0; h < K; h++) {
            H.body[i][V.body[i][h]] = 1;
        }
    }

    free_matrix(V);

    return H;
}


/**
    Collects data about information and check columns
*/
columns_metadata create_columns_metadata(int* check_set, int n, int k) {
	int R[n];
    int information_size = n;
    int i;

    for (i = 0; i < n; i++) {
        R[i] = 1;
    }

    for (i = 0; i < k; i++) {
        if (check_set[i] != -1) {
            R[check_set[i]] = 0;
            information_size--;
        }
    }

    int *information_set = (int *) malloc(information_size * sizeof(int));

    information_size = 0;

    for (i = 0; i < n; i++) {
        if (R[i] == 1) {
            information_set[information_size] = i;
            information_size++;
        }
    }

    check_set = (int*) realloc(check_set, (n - information_size) * sizeof(int));

    columns_metadata columns_mdata;
    columns_mdata.information_set = information_set;
    columns_mdata.check_set = check_set;
    columns_mdata.information_size = information_size;
    columns_mdata.check_size = n - information_size;

    return columns_mdata;
}

/**
    Gets max weight from weight matrix
*/
int get_max_weight(matrix weight_matrix) {
    int i = 0;
    int j = 0;
    int max = -1;
    for (i = 0; i < weight_matrix.rows; i++) {
        for (j = 0; j < weight_matrix.columns; j++) {
            if (max < weight_matrix.body[i][j]) {
                max = weight_matrix.body[i][j];
            }
        }
    }

    return max;
}

/**
    Gets array of pairs "weight-number"
*/
weight_number_pair* get_weight_number_pairs(matrix weight_matrix, int *num_of_weights_ptr) {
    int max_weight = -1;
    int i = 0;
    int j = 0;
    max_weight = get_max_weight(weight_matrix);
    int *weight_number_array = (int*)malloc(sizeof(int) * (max_weight + 1));

    for (i = 0; i < (max_weight + 1); i++) {
        weight_number_array[i] = 0;
    }

    for (i = 0; i < weight_matrix.rows; i++) {
        for (j = 0; j < weight_matrix.columns; j++) {
            weight_number_array[weight_matrix.body[i][j]]++;
        }
    }

    *num_of_weights_ptr = 0;
    for (i = 1; i < (max_weight + 1); i++) {
        if (weight_number_array[i] != 0) {
            *num_of_weights_ptr += 1;
        }
    }

   weight_number_pair *w_n_pair_array = (weight_number_pair*)malloc(sizeof(weight_number_pair) * *num_of_weights_ptr);
    int pair_index = 0;
    for (i = max_weight; i > 0; i--) {

        if (weight_number_array[i] != 0) {
            w_n_pair_array[pair_index].weight = i;
            w_n_pair_array[pair_index].number = weight_number_array[i];
            pair_index++;
        }
    }

    free(weight_number_array);

    return w_n_pair_array;
}

/**
    Generates polynomial using array of distances between ones
*/
int* generate_polynom_use_distances(int weight, int *distances, int submatrix_size) {

    int *polynom = NULL;
    polynom = (int*)calloc(submatrix_size, sizeof(int));
    int i = 0;

    for (i = 0; i < submatrix_size; i++) {
        polynom[i] = 0;
    }
    int index = 0;
    polynom[index] = 1;

    for (i = 0; i < weight - 1; i++) {
        index += distances[i];
        polynom[index] = 1;
    }

    return polynom;
}

/**
    Gets array of distances between ones
*/
int* get_distances_array(int *possible_distances_array, int num_of_distances, int weight, int submatrix_size) {

    int *distances_array = (int*)malloc(weight * sizeof(int));
    int i = 0;
    int j = 0;

    for (i = 0; i < weight; i++) {
        distances_array[i] = 0;
    }

    int distance_is_used = TRUE;
    while (distance_is_used) {
        distances_array[0] = possible_distances_array[(rand() % num_of_distances)];
        for (i = 1; i < weight - 1; i++) {
            int distance_not_unique = TRUE;
            while (distance_not_unique) {
                distance_not_unique = FALSE;
                distances_array[i] = possible_distances_array[(rand() % num_of_distances)];
                for (j = 0; j < i; j++) {
                    if (distances_array[i] == distances_array[j]) {
                        distance_not_unique = TRUE;
                        break;
                    }
                }
            }
        }

        int last = submatrix_size;
        for (i = 0; i < (weight - 1); i++) {
            last -= distances_array[i];
        }
        
        distance_is_used = FALSE;
        for (i = 0; i < weight - 1; i++) {
        	if (last == distances_array[i]) {
        		distance_is_used = TRUE;
        		break;
        	}
        }
        
        if (distance_is_used) continue;

		distance_is_used = TRUE;
        for (i = 0; i < num_of_distances; i++) {
            if (last == possible_distances_array[i]) {
                distance_is_used = FALSE;
                distances_array[weight - 1] = last;
                break;
            }
        }
    }
    return distances_array;
}

/**
    Deletes distances from array of possible distances
*/
int* delete_distances_from_array(int *possible_distances_array, int *num_of_possible_distances, int *distances_array, int weight) {
    int i = 0;
    int j = 0;
    int buf = 0;
    int last_index = *num_of_possible_distances - 1;
    for (i = 0; i < *num_of_possible_distances; i++) {
		for (j = 0; j < weight; j++) {
            if (possible_distances_array[i] == distances_array[j]) {
                buf = possible_distances_array[i];
                possible_distances_array[i] = possible_distances_array[last_index];
                possible_distances_array[last_index] = buf;
                last_index--;
                (*num_of_possible_distances)--;
            }
        }
    }

    possible_distances_array = (int*)realloc(possible_distances_array, sizeof(int) * (*num_of_possible_distances));
    return possible_distances_array;
}

/**
    Gets polynomial matrix
*/
void get_polynomial_matrix(matrix weight_matrix, int submatrix_size, int num_of_weights, weight_number_pair *w_n_pairs, int ***polynomial_matrix) {

    int *possible_distances_array = (int*)malloc(sizeof(int) * (submatrix_size - 1));
    int num_of_distances = submatrix_size - 1;
    int *num_of_distances_pointer = &num_of_distances;
    int i = 0;
    int j = 0;
    int k = 0;
    int l = 0;

    for (i = 0; i < num_of_distances; i++) {
        possible_distances_array[i] = i + 1;
    }

    for (i = 0; i < weight_matrix.rows; i++) {
        for (j = 0; j < weight_matrix.columns; j++) {
            for (k = 0; k < submatrix_size; k++) {
                polynomial_matrix[i][j][k] = 0;
            }
        }
    }

    int *current_distances_array;
    int *buf;
    for (i = 0; i < num_of_weights; i++) {
        k = 0;
        l = 0;
        for (j = 0; j < w_n_pairs[i].number; j++) {
            //MAGIC get distances array
            current_distances_array = get_distances_array(
                                        possible_distances_array,
                                        num_of_distances,
                                        w_n_pairs[i].weight,
                                        submatrix_size);
            //delete used distances from array
            possible_distances_array = delete_distances_from_array(
                                            possible_distances_array,
                                            num_of_distances_pointer,
                                            current_distances_array,
                                            w_n_pairs[i].weight);

            //searching for place for polynomial
            int polynom_generated = FALSE;
            for (; k < weight_matrix.rows; k++) {
                for (; l < weight_matrix.columns; l++) {
                    if (weight_matrix.body[k][l] == w_n_pairs[i].weight) {
                        //generate polynomial using distances
                        buf = NULL;
                        buf = generate_polynom_use_distances(
                                    w_n_pairs[i].weight,
                                    current_distances_array,
                                    submatrix_size);
                        int m = 0;
                        for (m = 0; m < submatrix_size; m++) {
                            polynomial_matrix[k][l][m] = buf[m];
                        }
                        free(buf);
                        buf = NULL;
                        polynom_generated = TRUE;
                        l++;
                        break;
                    }
            	}
                if (polynom_generated) {
                	if (l == weight_matrix.columns) {
                		k++;
                		l = 0;
					}
                    break;
                }
                if (l == weight_matrix.columns) l = 0;
            }
            free(current_distances_array);
        }
    }
    free(possible_distances_array);
}

/**
    Adds shifts to polynomial matrix
*/
void get_polynomial_matrix_with_shift(int ***polynomial_matrix, matrix weight_matrix, int submatrix_size) {

	int i = 0;
	int j = 0;
	int k = 0;
	int l = 0;

	for (i = 1; i < weight_matrix.rows; i++) {
		for (j = 0; j < weight_matrix.columns; j++) {
			int shift = 0;
			int position_exists = TRUE;

			while (position_exists) {
				shift++;
				position_exists = FALSE;
				for (k = 0; k < i; k++) {
					for (l = 0; l < submatrix_size; l++) {
						if ((polynomial_matrix[i][j][l] == 1) && (polynomial_matrix[k][j][(l + shift) % submatrix_size] == 1)) {

							position_exists = TRUE;
							break;
						}
					}
					if (position_exists) {
						break;
					}
				}
			}

			int polynom_copy[submatrix_size];

			for (l = 0; l < submatrix_size; l++) {
				polynom_copy[(l + shift) % submatrix_size] = polynomial_matrix[i][j][l];
			}

			for (l = 0; l < submatrix_size; l++) {
				polynomial_matrix[i][j][l] = polynom_copy[l];
			}
		}
	}
}

/**
    Creates H matrix using polynomial matrix with shifts
*/
matrix create_H_matrix_use_polynomial_matrix_with_shifts(int ***polynomial_matrix, matrix weight_matrix, int submatrix_size) {

	matrix H_matrix = create_empty_matrix(weight_matrix.rows * submatrix_size, weight_matrix.columns * submatrix_size);
	int i = 0;
	int j = 0;
	int k = 0;
	int l = 0;

	for (i = 0; i < weight_matrix.rows; i++) {
		for (j = 0; j < weight_matrix.columns; j++) {
			int *polynomial = polynomial_matrix[i][j];
			for (k = 0; k < submatrix_size; k++) {
				for (l = 0; l < submatrix_size; l++) {
					H_matrix.body[i * submatrix_size + k][j * submatrix_size + l] = polynomial[((l - k) + submatrix_size) % submatrix_size];
				}
			}
		}
	}

	return H_matrix;
}

/**
    Creates H matrix of NEW_CODE type
*/
matrix create_H_matrix_of_new_code(matrix weight_matrix, int submatrix_size) {

    int ***polynomial_matrix = create_three_dimensional_array(weight_matrix.rows, weight_matrix.columns, submatrix_size);
    int num_of_weights = 0;
    weight_number_pair *w_n_pairs = get_weight_number_pairs(weight_matrix, &num_of_weights);
    get_polynomial_matrix(weight_matrix, submatrix_size, num_of_weights, w_n_pairs, polynomial_matrix);
	get_polynomial_matrix_with_shift(polynomial_matrix, weight_matrix, submatrix_size);
    matrix H = create_H_matrix_use_polynomial_matrix_with_shifts(polynomial_matrix, weight_matrix, submatrix_size);
    free_three_dimensional_array(polynomial_matrix, weight_matrix.rows, weight_matrix.columns);

    return H;
}


/**
    Creates LDPC code structure
*/
ldpc create_ldpc(code_type type, int J, int K, int M, matrix weight_matrix) {
	matrix H;
    if (type != NEW_CODE) {
        H = create_H_rand(type, J, K, M);
    } else {
    	H = create_H_matrix_of_new_code(weight_matrix, M);
	}

	matrix H_copy = copy_matrix(H);
    int* check_set = perform_gauss_elimination(H_copy);
    columns_metadata columns_mdata = create_columns_metadata(check_set, K * M, J * M);

    matrix cutted_H;
    switch (type)
    {
        case GALLAGER:
            cutted_H = copy_matrix_part(H_copy, M * J - J + 1, K * M);
            break;
        case RU_CODE:
            cutted_H = copy_matrix(H_copy);
            break;
        case NEW_CODE:
        	cutted_H = copy_matrix(H_copy);
        default:
            cutted_H = create_random_matrix(J * M, K * M);
    }

    matrix G = create_G_from_H_matrix(cutted_H, columns_mdata);

    free_matrix(cutted_H);

    ldpc ldpc_object;
    ldpc_object.G = G;
    ldpc_object.H = H;
    ldpc_object.columns_mdata = columns_mdata;
    ldpc_object.n = G.columns;
    ldpc_object.k = G.rows;
    ldpc_object.r = H.rows;
    ldpc_object.systematic_r = cutted_H.rows;
    ldpc_object.systematic_H = H_copy;
    ldpc_object.C = get_non_zero_column_data(H);
    ldpc_object.V = get_non_zero_column_data(transpose_matrix(H));

    return ldpc_object;
}

/**
    Destroys LDPC code structure
*/
void free_ldpc(ldpc ldpc_object) {
    free_matrix(ldpc_object.G);
    free_matrix(ldpc_object.H);
	free(ldpc_object.columns_mdata.check_set);
	free(ldpc_object.columns_mdata.information_set);
}

/**
    Gets data about indices of nonero elements of H matrix
*/
indices_of_nonzero_elements get_non_zero_column_data(matrix matrix_object) {
    indices_of_nonzero_elements non_zero_column_data;
    non_zero_column_data.element_data = (int **)malloc(matrix_object.columns * sizeof(int *));
    non_zero_column_data.element_length = (int *)malloc(matrix_object.columns * sizeof(int));
    int buf_column[matrix_object.rows];
    int i, j;

    for (i = 0; i < matrix_object.columns; i++) {
        int num_ones = 0;
        for (j = 0; j < matrix_object.rows; j++) {
            if (matrix_object.body[j][i] != 0) {
                buf_column[num_ones] = j;
                num_ones++;
            }
        }

        non_zero_column_data.element_length[i] = num_ones;
        non_zero_column_data.element_data[i] = (int *)malloc(num_ones * sizeof(int));

        for (j = 0; j < num_ones; j++) {
            non_zero_column_data.element_data[i][j] = buf_column[j];
        }
    }

    return non_zero_column_data;
}
