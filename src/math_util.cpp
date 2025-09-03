#include <stdio.h>
#include "math_util.h"


void assign_value(matrix_t *m, int row, int col, float value) {
    int dimCol, idx;
    dimCol = m->dimCol;
    idx = row * dimCol + col;
    m->data[idx] = value;
}


float get_value(matrix_t *m, int row, int col)
{
    int dimCol, idx;
    dimCol = m->dimCol;
    idx = row * dimCol + col;
    return m->data[idx];
}

/**
 * Multiply left and right matrices
 * Dimensions must match.

 * Mismatching dimensions leaves all matrices as is and sets and error code in
 * the errorcode argument
 * It is the responsibility of the caller to check this
 * 
 */
int matmul(matrix_t *left, matrix_t *right, matrix_t *result, int *errorcode)
{

    if (!(left->dimCol == right->dimRow))
    {
        fprintf(stderr, "matrix multiplications mismatch with shapes (%dx%d) * (%dx%d), errorcode %d\n",
        left->dimRow, left->dimCol, right->dimRow, right->dimCol, MATMUL_DIMENSION_MISMATCH_ERROR);
        *errorcode = MATMUL_DIMENSION_MISMATCH_ERROR;
        return 0;
    }

    float res;
    for (int row = 0; row < result->dimRow; row++)
    {
        for (int col = 0; col < result->dimCol; col++)
        {
            res = 0.0f;
            for (int i = 0; i < left->dimRow; i++)
            {
                res += get_value(left, row, i) * get_value(right, i, col);
            }
            assign_value(result, row, col, res);
        }
    }
    return 1;
}