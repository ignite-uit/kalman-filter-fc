#include <stdio.h>
#include "math_util.h"

// https://www.andreinc.net/2021/01/20/writing-your-own-linear-algebra-matrix-library-in-c#retrieving--selecting-a-column
void set_val(matrix_t *m, int row, int col, float value)
{
    int numCol, idx;
    numCol = m->numCol;
    idx = row * numCol + col;
    m->data[idx] = value;
}

float get_value(matrix_t *m, int row, int col)
{
    int numCol, idx;
    numCol = m->numCol;
    idx = row * numCol + col;
    return m->data[idx];
}

void pprint_matrix(matrix_t *A) {
    int N = A->numRow;
    int M = A->numCol;
    int i=0, j=0;
    for (i=0; i < N; i++) {
        printf("|");
        for (j=0; j < M; j++) {
            printf(" %8.4f", get_value(A, i, j));
        }
        printf(" |\n");
    }
}

int copy_matrix(matrix_t *matA, matrix_t *matB) {
    int N = matA->numRow, M = matA->numCol;

    if (!(N == matB->numRow && M == matB->numCol)) {
        fprintf(stderr, "attemt to copy matrix with mismatching dimensions in target matrix");
        return 0;
    }

    float *A = matA->data;
    float *B = matB->data;
    for (int j = 0; j < N * M; j++)
    {
        B[j] = A[j];
    }
    return 1;
}

void clear_matrix(matrix_t *mat) {
    int N = mat->numRow, M = mat->numCol;
    float *m = mat->data;
    for (int j = 0; j < N * M; j++)
    {
        m[j] = 0.0f;
    }
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

    if (!(left->numCol == right->numRow))
    {
        fprintf(stderr, "matrix multiplications mismatch with shapes (%dx%d) * (%dx%d), errorcode %d\n",
                left->numRow, left->numCol, right->numRow, right->numCol, MATMUL_DIMENSION_MISMATCH_ERROR);
        *errorcode = MATMUL_DIMENSION_MISMATCH_ERROR;
        return 0;
    }

    float res;
    for (int row = 0; row < result->numRow; row++)
    {
        for (int col = 0; col < result->numCol; col++)
        {
            res = 0.0f;
            for (int i = 0; i < left->numRow; i++)
            {
                res += get_value(left, row, i) * get_value(right, i, col);
            }
            set_val(result, row, col, res);
        }
    }
    return 1;
}

int matadd(matrix_t *a, matrix_t *b, matrix_t *result, int *errorcode)
{
    int numCol = result->numCol;
    int numRow = result->numRow;

    if (!((a->numCol == b->numCol && b->numCol == numCol) &&
          (a->numRow == b->numRow && b->numRow == numRow)))
    {
        fprintf(stderr, "matrix addition shape mismatch with shapes"
                        "(%dx%d) * (%dx%d) = (%dx%d), errorcode %d\n",
                a->numRow, a->numCol, b->numRow, b->numCol, result->numRow, result->numCol, MATADD_DIMENSION_MISMATCH_ERROR);
        *errorcode = MATADD_DIMENSION_MISMATCH_ERROR;
        return 0;
    }

    int row, col;
    float a_rc, b_rc;
    for (row = 0; row < numRow; row++)
    {
        for (col = 0; col < numCol; col++)
        {
            a_rc = get_value(a, row, col);
            b_rc = get_value(b, row, col);
            set_val(result, row, col, a_rc + b_rc);
        }
    }
    return 1;
}

int matsub(matrix_t *a, matrix_t *b, matrix_t *result, int *errorcode)
{
    int numCol = result->numCol;
    int numRow = result->numRow;

    if (!((a->numCol == b->numCol && b->numCol == numCol) &&
          (a->numRow == b->numRow && b->numRow == numRow)))
    {
        fprintf(stderr, "matrix subtraction shape mismatch with shapes"
                        "(%dx%d) * (%dx%d) = (%dx%d), errorcode %d\n",
                a->numRow, a->numCol, b->numRow, b->numCol, result->numRow, result->numCol, MATADD_DIMENSION_MISMATCH_ERROR);
        *errorcode = MATADD_DIMENSION_MISMATCH_ERROR;
        return 0;
    }

    int row, col;
    float a_rc, b_rc;
    for (row = 0; row < numRow; row++)
    {
        for (col = 0; col < numCol; col++)
        {
            a_rc = get_value(a, row, col);
            b_rc = get_value(b, row, col);
            set_val(result, row, col, a_rc - b_rc);
        }
    }
    return 1;
}

/**
 * swap rows i and j in place
 */
void swaprows_inplace(matrix_t *m, int i, int j)
{
    float elem_i, elem_j;
    for (int col = 0; col < m->numCol; col++)
    {
        elem_i = get_value(m, i, col);
        elem_j = get_value(m, j, col);
        set_val(m, i, col, elem_j);
        set_val(m, j, col, elem_i);
    }
}

/**
 * swap columns i and j in place
 */
void swapcols_inplace(matrix_t *m, int i, int j)
{
    float elem_i, elem_j;
    for (int row = 0; row < m->numRow; row++)
    {
        elem_i = get_value(m, row, i);
        elem_j = get_value(m, row, j);
        set_val(m, row, i, elem_j);
        set_val(m, row, j, elem_i);
    }
}

void mult_mat_scal(matrix_t *mat, float scalar)
{
    int N = mat->numRow, M = mat->numCol;
    float *m = mat->data;
    for (int j = 0; j < N * M; j++)
    {
        m[j] = m[j] * scalar;
    }
}

/**
 * take inverse of 3x3 matrix A and store it in matrix invA
 *
 */
int inv3x3(matrix_t *A, matrix_t *invA)
{
    if (!(A->numRow == 3 && A->numCol==3 && invA->numCol==3 && invA->numRow == 3)) {
        fprintf(stderr, "shapes of matrix A and invA must both be 3x3\n");
    }

    float a11, a12, a13, a21, a22, a23, a31, a32, a33; // matrix elements
    float m11, m12, m13, m21, m22, m23, m31, m32, m33; // minors
    float c11, c12, c13, c21, c22, c23, c31, c32, c33; // cofactors
    float det3x3;
    // clang-format off
    a11 = get_value(A, 0, 0); a12 = get_value(A, 0, 1); a13 = get_value(A, 0, 2);
    a21 = get_value(A, 1, 0); a22 = get_value(A, 1, 1); a23 = get_value(A, 1, 2);
    a31 = get_value(A, 2, 0); a32 = get_value(A, 2, 1); a33 = get_value(A, 2, 2);

    // minors
    m11 = a22*a33 - a23*a23; m12 = a21*a33 - a31*a23; m13 = a21*a32 - a31*a22;
    m21 = a12*a33 - a32*a13; m22 = a11*a33 - a31*a13; m23 = a11*a32 - a31*a12; 
    m31 = a12*a23 - a22*a13; m32 = a11*a23 - a21*a13; m33 = a11*a22 - a21*a12;


    // cofactor c_ij = (-1)^(i + j)*m_ij
    c11 =      m11, c12 = (-1)*m12, c13 =      m13;
    c21 = (-1)*m21, c22 =      m22, c23 = (-1)*m23;
    c31 =      m31, c32 = (-1)*m32, c33 =      m33;

    det3x3 = a11 * c11 + a12 * c12 + a13 * c13;
    printf("detA = %f\n", det3x3);

    if (det3x3 < 1E-5)
    {
        // noninvertible matrix!
        fprintf(stderr, "Error: taking inverse of non-invertible matrix!");
        return -1;
    }

    // form the transposed cofactor matrix C^T: index i,j --> j,i
    set_val(invA, 0, 0, c11); set_val(invA, 0, 1, c21); set_val(invA, 0, 2, c31);
    set_val(invA, 1, 0, c12); set_val(invA, 1, 1, c22); set_val(invA, 1, 2, c32);
    set_val(invA, 2, 0, c13); set_val(invA, 2, 1, c23); set_val(invA, 2, 2, c33);
    // clang-format on

    mult_mat_scal(invA, 1 / det3x3);
    return 0;
}