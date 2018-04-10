#include "matrix.h"

#ifndef LDPC_TESTER
#define LDPC_TESTER

typedef struct{
	float min, max, step;
} SNR_interval;

matrix create_random_message(int length);

#endif
