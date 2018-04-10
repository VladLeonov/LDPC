#include "ldpc_tester.h"

#include <stdlib.h>

matrix create_random_message(int length) {
    int i;
    matrix message = create_zero_matrix(1, length);
    for (i = 0; i < length; i++) {
        message.body[0][i] = rand() % 2;
    }
    return message;
}
