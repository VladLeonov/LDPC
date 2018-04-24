#include <math.h>
#include "subsidary_math.h"
#include "stdlib.h"

void fill_with_permutation(int *x, int n) {
	int i, j, temp;
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

int get_indexes_of_common_elements(int *arr_a, int *arr_b, int *result, int len_a, int len_b) {
    int i, j, result_length;
    j = 0;
    result_length = 0;

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


float log_exp(float x) {
    x = fmaxf(fminf(x, 19.07), 0.01);
    return log((exp(x) - 1) / (exp(x) + 1));
}


float log_tahn(float value) {
    float t = exp(fabs(value));
    return log((t + 1)/(t - 1));
}


int sign(float value) {
    if (value > 0) {
        return 1;
    } else if (value < 0) {
        return -1;
    } else {
        return 0;
    }
}


float randn() {
	float U1, U2, W, mult;
	static float X1, X2;
	static int call = 0;

	if (call == 1) {
		call = !call;
		return X2;
	}

	do {
		U1 = -1 + ((float) rand () / RAND_MAX) * 2;
		U2 = -1 + ((float) rand () / RAND_MAX) * 2;
		W = pow (U1, 2) + pow (U2, 2);
	} while ((W >= 1) || (W == 0));

	mult = sqrt ((-2 * log (W)) / W);
	X1 = U1 * mult;
	X2 = U2 * mult;

	call = !call;

	return X1;
}
