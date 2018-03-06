#include "matrix.h"
#include "ldpc.h"

int* gauss_elimination(matrix G) {
	int k = G.rows;
	int n = G.columns;
    int* information_set = (int*) malloc(k * sizeof(int));
    int i, j;
    for (i = 0; i < k; i++){
        information_set[i] = -1;
    }
    
    i = 0;
    int column_number = 0;

    while (i < k) {
        
        //throw away all-zero row
        while (sum_rows(G, i) == 0) {
            for (j = i; j < (k - 1); j++) {
                G.body[j] = G.body[j + 1];
            }
            k--;
            free(G.body[k]);
            G.rows--;

            information_set[k] = -1;
            
            if (i >= k) {
                break;
            }
        }
        
        if (i >= k) {
            break;
        }
        
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
            char *buff = G.body[i];
            G.body[i] = G.body[j];
            G.body[j] = buff;

            i++;
        }

        column_number++;  // keep the same row, another column
    }
    return information_set;
}

ldpc create_systematic_view(matrix G_old) {

    matrix G = copy_matrix(G_old);

    int k = G.rows;
    int n = G.columns;
    
    int* information_set = gauss_elimination(G);

    // transposition
    char R[n];
    int check_size = n;
    int i, j;
    for (i = 0; i < n; i++) {
        R[i] = 1;
    }
    for (i = 0; i < k; i++) {
        if (information_set[i] != -1) {
            R[information_set[i]] = 0;
            check_size--;
        }
    }
    
    int *check_set = (int *) malloc(check_size);
    check_size = 0;
    for (i = 0; i < n; i++) {
        if (R[i] == 1) {
            check_set[check_size] = i;
            check_size++;
        }
    }

    matrix H = create_zero_matrix(n-k, n);
    matrix P = create_empty_matrix(G.rows, check_size);
    for (j = 0; j < check_size; j++) {
        for (i = 0; i < G.rows; i++) {
            P.body[i][j] = G.body[i][check_set[j]];
        }
    }
    
    matrix PT = transpose_matrix(P);
    free_matrix(P);
    for (j = 0; j < n - check_size; j++) {
        for (i = 0; i < H.rows; i++) {
            H.body[i][information_set[j]] = PT.body[i][j];
        }
    }
    
    matrix U = create_unit_matrix(n - k);
    for (j = 0; j < check_size; j++) {
        for (i = 0; i < H.rows; i++) {
            H.body[i][check_set[j]] = U.body[i][j];
        }
    }
    
    ldpc ldpc_object;
    ldpc_object.G = G;
    ldpc_object.H = H;
    ldpc_object.check_set = check_set;
    ldpc_object.information_set = information_set;
    ldpc_object.k = G.rows;
    ldpc_object.n = G.columns;
    ldpc_object.check_size = check_size;
    ldpc_object.information_size = n - check_size;

    return ldpc_object;
}

int sum_rows(matrix G, int row_index) {
    int i;
    int sum = 0;
    for (i = 0; i < G.columns; i++) {
        sum += G.body[row_index][i];
    }
    return sum;
}

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
