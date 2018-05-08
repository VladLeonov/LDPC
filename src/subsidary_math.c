/**
    LDPC
    subsidary_math.c
    Purpose: Ñontains subsidary mathematical operations 
	for using by other methods.

    @author Leonov V.R.
    @version 24.04.18
*/


#include <math.h>
#include <stdlib.h>

#include "subsidary_math.h"


#define TRUE !0
#define FALSE 0


//Fills int array with random permutation.
void fill_with_permutation(int *x, int n) {
	int i = 0, j = 0;
	int temp = 0;
	
    for (i = 0; i < n; i++) {
        x[i] = i;
    }
    
    for (i = n-1; i > 0; i--) {
        j = rand() % (i+1);
        temp = x[j];
        x[j] = x[i];
        x[i] = temp;
    }
}


//Finds common elements of two arrays 
//and returns their indexes in the first array.
int get_indexes_of_common_elements(int *arr_a, int *arr_b, int *result,
                                   int len_a, int len_b) {
    int i = 0, j = 0;
	int result_length = 0;

    for (i = 0; i < len_a; i++) {
        for (j = 0; j < len_b; j++) {
            if (arr_a[i] == arr_b[j]) {
                result[result_length] = i;
                result_length++;
                break;
            } else if (arr_a[i] < arr_b[j]) {
                break;
            } else {
                continue;
            }
        }
    }

    return result_length;
}


//Calculate log((exp(x) - 1) / (exp(x) + 1)).
float log_exp(float x) {
    x = fmaxf(fminf(x, 19.07), 0.01);
    return log((exp(x) - 1) / (exp(x) + 1));
}


//Calculate log((exp(x) + 1) / (exp(x) - 1)).
float log_tahn(float value) {
    float t = exp(fabs(value));
    return log((t + 1)/(t - 1));
}


//Determines the sign of a value.
int sign(float value) {
    if (value > 0) {
        return 1;
    } else if (value < 0) {
        return -1;
    } else {
        return 0;
    }
}


//Generates a random number from the normal distribution.
float randn() {
	static float X1 = 0, X2 = 0;
	static char call = FALSE;

	if (call == TRUE) {
		call = !call;
		return X2;
	}

    float U1 = 0, U2 = 0;
    float W = 0;

	do {
		U1 = -1 + ((float) rand () / RAND_MAX) * 2;
		U2 = -1 + ((float) rand () / RAND_MAX) * 2;
		W = pow (U1, 2) + pow (U2, 2);
	} while ((W >= 1) || (W == 0));
	
	float mult = 0;

	mult = sqrt ((-2 * log (W)) / W);
	X1 = U1 * mult;
	X2 = U2 * mult;

	call = !call;

	return X1;
}


//Creates three-dimensional array filled with zeros.
int*** create_three_dimensional_array(int x, int y, int z) {
	int ***array = (int***) malloc(x * sizeof(int**));
	int i, j, k;
    for (i = 0; i < x; i++) {
        array[i] = (int**) malloc(y * sizeof(int*));
        for (j = 0; j < y; j++) {
            array[i][j] = (int*) malloc(z * sizeof(int));
            for (k = 0; k < z; k++) {
                array[i][j][k] = 0;
            }
        }
    }
    return array;
}


//Frees the memory occupied by the three-dimensional array.
void free_three_dimensional_array(int*** array, int x, int y) {
	int i, j, k;
    for (i = 0; i < x; i++) {
        for (j = 0; j < y; j++) {
            free(array[i][j]);
        }
        free(array[i]);
    }
    free(array);
}
