#include "ldpc.h"
#include "matrix.h"
#include "math.h"

int main() {

    /*matrix message = create_empty_matrix(1, 20);
    int i;
    for (i = 0; i < 20; i++) {
        message.body[0][i] = 'a' + i;
    }

    ldpc ldpc_object = create_ldpc(7, 4);
    matrix encoded = encode(ldpc_object, message);
    matrix decoded = decode(ldpc_object, encoded);

    for (i = 0; i < 20; i++) {
        printf("%c ", message.body[0][i]);
    }*/

    int array[3][5] = {{0, 0, 1, 1, 1},
                       {1, 0, 1, 1, 0},
                       {1, 0, 0, 1, 1}};
    matrix G = array_to_matrix(3, 5, array);
    ldpc ldpc_object = create_systematic_view(G);

    printf("Not systematic G = \n");
    print_matrix(G); 
    printf("\n");
    print_ldpc(ldpc_object);
    
    printf("G * HT = \n");
    print_matrix(multiply_matrices(ldpc_object.G, transpose_matrix(ldpc_object.H)));
    printf("\n");

    system("pause");
    return 0;
}
