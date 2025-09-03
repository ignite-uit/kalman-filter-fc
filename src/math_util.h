
#define MATMUL_DIMENSION_MISMATCH_ERROR 1

typedef struct matrix
{
    int dimCol;
    int dimRow;
    float *data;
} matrix_t;

typedef struct vector
{
    int dim;
    float *data;
} vector_t;

void assign_value(matrix_t *m, int row, int col, float value);

float get_value(matrix_t *m, int row, int col);

/**
 * Multiply left and right matrices
 * Dimensions must match.

 * Mismatching dimensions leaves all matrices as is and sets and error code in
 * the errorcode argument
 * It is the responsibility of the caller to check this
 * 
 */
int matmul(matrix_t *left, matrix_t *right, matrix_t *result, int *errorcode);