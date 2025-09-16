#include <stdio.h>
#include "math_util.h"


// static memory allocations for matrices
float Id_data[3 * 3] = {0.0f};
float A_data[3 * 3] = {0.0f};
float D_data[3 * 3] = {0.0f};
float invA_data[3 * 3] = {0.0f};
float T_data[3 * 3] = {0.0f};
float invT_data[3 * 3] = {0.0f};

float S_data[3 * 3] = {0.0f}; // scratch matrix
float S2_data[3 * 3] = {0.0f}; // scratch matrix

matrix_t A, D, invA, Id, T, invT, S, S2;


void initialize()
{
    Id.data = Id_data;
    A.data = A_data;
    D.data = D_data;
    invA.data = invA_data;
    T.data = T_data;
    invT.data = invT_data;
    S.data = S_data;
    S2.data = S2_data;

    Id.numRow = 3;
    Id.numCol = 3;
    set_val(&Id, 0, 0, 1);
    set_val(&Id, 1, 1, 1);
    set_val(&Id, 2, 2, 1);

    // 
    D.numRow = 3;
    D.numCol = 3;
    set_val(&D, 0, 0, 5);
    set_val(&D, 1, 1, 5);
    set_val(&D, 2, 2, 5);

    S.numRow = 3;
    S.numCol = 3;
    S2.numRow = 3;
    S2.numCol = 3;


    // set some linearly independent columns
    // test case: A != inv(A) and A * inv(A) = Id
    A.numRow = 3;
    A.numCol = 3;
    set_val(&A, 0, 3, 1);
    set_val(&A, 1, 2, 1);
    set_val(&A, 2, 1, 1);

    set_val(&A, 0, 1, 1);
    set_val(&A, 1, 1, 2);
    set_val(&A, 2, 1, 3);

    set_val(&A, 0, 2, 4);
    set_val(&A, 1, 2, 0);
    set_val(&A, 2, 2, -9);
    
}

int main()
{
    initialize();

    pprint_matrix(&Id);
    printf("\n");
    printf("Matrix D\n");
    pprint_matrix(&D);

    int e;
    inv3x3(&D, &S, &e);
    printf("inverse D\n");
    pprint_matrix(&S);

    printf("matmul D * inv(D)\n");
    matmul(&D, &S, &S2, &e);
    pprint_matrix(&S2);


    clear_matrix(&S);
    clear_matrix(&S2);
    printf("\n");
    printf("\n");
    printf("matrix A\n");
    pprint_matrix(&A);
    printf("\n");

    inv3x3(&A, &S, &e);
    printf("inverse A\n");
    pprint_matrix(&S);

    printf("\n");
    matmul(&A, &S, &S2, &e);
    printf("\n");
    printf("A * inv(A): \n");
    pprint_matrix(&S2);






}