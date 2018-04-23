#include <math.h>
#include <stdlib.h>
#include "decoder.h"
#include "matrix.h"
#include "subsidary_math.h"

#define TRUE !0
#define FALSE 0
#define MAXITER 50

float* map_sp(float y[], int length) {
    // LLR domain
    int hard[length], synd = 0;
    int i;
    for (i = 0; i < length; i++) {
        if (y[i] < 0) {
            hard[i] = 1;
            synd++;
        } else {
            hard[i] = 0;
        }
    }
    synd %= 2;

    for (i = 0; i < length; i++) {
        hard[i] = (hard[i] + synd) % 2;
    }

    float alogpy[length], sum_alogpy = 0;
    for (i = 0; i < length; i++) {
        alogpy[i] = log_exp(fabs(y[i]));
        sum_alogpy += alogpy[i];
    }

    float *soft_out = (float*) malloc(length * sizeof(float));
    for (i = 0; i < length; i++) {
    	soft_out[i] = (2 * hard[i] - 1) * log_exp(alogpy[i] - sum_alogpy);
	}
    return soft_out;
}

char check_syndrome(matrix hard, int r, indices_of_nonzero_elements V) {
	int i, j, syn;
    for (i = 0; i < r; i++) {
    	syn = 0;
    	for (j = 0; j < V.element_length[i]; j++) {
    		syn += hard.body[0][V.element_data[i][j]];
		}
        syn %= 2;
        if (syn > 0) {
            return FALSE;
        }
    }
    return TRUE;
}

matrix get_hard_from_soft(float soft[], int length) {
    matrix result = create_zero_matrix(1, length);

    int i;
    for (i = 0; i < length; i++) {
        if (soft[i] < 0) {
            result.body[0][i] = 1;
        } else {
            result.body[0][i] = 0;
        }
    }

    return result;
}

int flooding(ldpc ldpc_object, float *soft, matrix *hard_solution) {
	int i, j;
    int n = ldpc_object.H.columns;
    int r = ldpc_object.H.rows; // number of parity checks
    indices_of_nonzero_elements V = ldpc_object.V;
    indices_of_nonzero_elements C = ldpc_object.C;

    int rw[r];
    for (i = 0; i < r; i++) {
    	rw[i] = V.element_length[i];
    }
    int cw[n];
    for (i = 0; i < n; i++) {
    	cw[i] = C.element_length[i];
    }

    float Z[r][n];  // current LLRS
    for (i = 0; i < r; i++) {
    	for (j = 0; j < n; j++) {
	    	Z[i][j] = 0;
	    }
    }

    float soft_out[n];
	for (i = 0; i < n; i++) {
        soft_out[i] = soft[i];
    }

    matrix hard = get_hard_from_soft(soft, n);

    if (check_syndrome(hard, r, V) == TRUE) {
    	hard_solution->body = hard.body;
    	hard_solution->columns = hard.columns;
    	hard_solution->rows = hard.rows;
        return 0;
    }

    for (i = 0; i < r; i++) {
    	for (j = 0; j < rw[i]; j++) {
    		Z[i][V.element_data[i][j]] = soft[V.element_data[i][j]];
    	}
    }

	int steps;
    for (steps = 1; steps <= MAXITER; steps++) {
        // loop over checks
        float* soft_buffer;
        float y[n];
        for (i = 0; i < r; i++) {
        	for (j = 0; j < rw[i]; j++) {
	    		y[j] = Z[i][V.element_data[i][j]];
	    	}
	    	soft_buffer = map_sp(y, rw[i]);

	        for (j = 0; j < rw[i]; j++) {
	    		Z[i][V.element_data[i][j]] = soft_buffer[j];
	    	}
	    	free(soft_buffer);
        }

        // symbol nodes
        for (i = 0; i < n; i++) {
            // prob domain
            float sum_Z = 0;
            for (j = 0; j < cw[i]; j++) {
            	sum_Z += Z[C.element_data[i][j]][i];
            }
            soft_out[i] = soft[i] + sum_Z;
            for (j = 0; j < cw[i]; j++) {
				Z[C.element_data[i][j]][i] = soft_out[i] - Z[C.element_data[i][j]][i];
			}
        }

        hard = get_hard_from_soft(soft_out, n);

        if (check_syndrome(hard, r, V) == TRUE) {
        	hard_solution->body = hard.body;
    		hard_solution->columns = hard.columns;
    		hard_solution->rows = hard.rows;
            return steps;
        }
    }

    hard_solution->body = hard.body;
    hard_solution->columns = hard.columns;
    hard_solution->rows = hard.rows;

    return -MAXITER;
}

int decode_belief_propogandation(ldpc ldpc_object, float *y, matrix *hard_solution, char use_non_zero_data) {
    matrix H = ldpc_object.H;
    int r = H.rows;
    int n = H.columns;
    float Z[r][n];
    float L[r][n];
    indices_of_nonzero_elements C = ldpc_object.C;
    indices_of_nonzero_elements V = ldpc_object.V;

    int i, j;
    for (i = 0; i < r; i++) {
        for (j = 0; j < n; j++) {
            Z[i][j] = 0;
            L[i][j] = 0;
        }
    }

    float soft[n];
    for (i = 0; i < n; i++) {
        soft[i] = y[i];
    }

    matrix hard = get_hard_from_soft(soft, n);

    matrix syndrome = count_syndrome(ldpc_object, hard, use_non_zero_data);
    int iter = 0;
    float L_element = 0;
    int index_C = 0;
    int k;

    while ((iter < MAXITER) && (sum_row_elements(syndrome, 0) != 0)) {
        //H columns processing
        float sum_Z;
        for (i = 0; i < n; i++) {
        	sum_Z = 0;
            for (j = 0; j < C.element_length[i]; j++) {
                sum_Z += Z[C.element_data[i][j]][i];
            }
            for (j = 0; j < C.element_length[i]; j++) {
                index_C = C.element_data[i][j];
                L[index_C][i] = y[i] + sum_Z - Z[index_C][i];
            }
        }

        //H rows processing
		float a, b;
		int index_C;
		for (i = 0; i < n; i++) {
		    for (index_C = 0; index_C < C.element_length[i]; index_C++) {
		        int j = C.element_data[i][index_C];
		        a = 1;
		        b = 0;
		        for (k = 0; k < V.element_length[j]; k++) {
		        	if (V.element_data[j][k] != i) {
			        	L_element = L[j][V.element_data[j][k]];
			            a *= sign(L_element);
			            b += log_tahn(L_element);
					}
		        }
		        Z[j][i] = fmaxf(fminf(a * log_tahn(b), 19.07), -19.07);
		    }
		}

        //result forming
        for (i = 0; i < n; i++) {
            soft[i] = y[i] + sum_coloumn_elements(n, Z, i, r);
        }

        free_matrix(hard);
        hard = get_hard_from_soft(soft, n);
        free_matrix(syndrome);
        syndrome = count_syndrome(ldpc_object, hard, use_non_zero_data);
        iter++;
    }

    hard_solution->body = hard.body;
    hard_solution->columns = hard.columns;
    hard_solution->rows = hard.rows;

    free_matrix(syndrome);

    return iter;
}

matrix count_syndrome(ldpc ldpc_object, matrix codedword, char use_non_zero_data) {

    matrix syndrome;
    if (use_non_zero_data == FALSE) {

    	matrix H_transposed = transpose_matrix(ldpc_object.H);
	    syndrome = multiply_matrices(codedword, H_transposed);
	    free_matrix(H_transposed);

    } else {

    	syndrome = create_zero_matrix(1, ldpc_object.r);
    	int i, j;
    	for (i = 0; i < ldpc_object.r; i++) {
    		for (j = 0; j < ldpc_object.V.element_length[i]; j++) {
    			syndrome.body[0][i] ^= codedword.body[0][ldpc_object.V.element_data[i][j]];
			}
		}

	}

    return syndrome;
}
