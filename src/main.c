#include "ldpc.h"
#include "matrix.h"
#include "math.h"
#include "ldpc_tester.h"

#include <stdlib.h>
#include <stdio.h>

#define TRUE !0
#define FALSE 0

int main() {

    int J = 9, K = 18, M = 8;
    ldpc ldpc_object = create_ldpc(RU_code, J, K, M);
    
    SNR_interval SNR = {0.0, 1.0, 0.1};
    FILE *file = fopen("decoding_simulation.txt", "w");
    decoding_simulation(ldpc_object, SNR, file);
    fclose(file);

    //int message_size = M * (K - J) + J - 1;
    
    system("pause");

    free_ldpc(ldpc_object);

    return 0;
}
