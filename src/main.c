#include "ldpc.h"
#include "matrix.h"
#include "math.h"
#include "ldpc_tester.h"

#include <stdlib.h>
#include <stdio.h>

#define TRUE !0
#define FALSE 0

int main() {

    int J = 3, K = 6, M = 3;//4, 8, 32
    ldpc ldpc_object = create_ldpc(RU_code, J, K, M);

    SNR_interval SNR = {3.8, 3.9, 0.005};//{3.8, 3.9, 0.005};
    FILE *file = fopen("decoding_simulation.txt", "w");
    decoding_simulation(ldpc_object, SNR, file);
    fclose(file);

    //int message_size = M * (K - J) + J - 1;

    system("pause");
    free_ldpc(ldpc_object);

    return 0;
}
