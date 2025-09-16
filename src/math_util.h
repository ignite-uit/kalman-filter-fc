
#define MATMUL_DIMENSION_MISMATCH_ERROR 1
#define MATADD_DIMENSION_MISMATCH_ERROR 2

typedef struct matrix
{
    int numCol;
    int numRow;
    float *data;
} matrix_t;

typedef struct vector
{
    int dim;
    float *data;
} vector_t;

/** Pretty print a matrix A */
void pprint_matrix(matrix_t *A);

void set_val(matrix_t *m, int row, int col, float value);

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



int copy_matrix(matrix_t *matA, matrix_t *matB);
void clear_matrix(matrix_t *mat);


/**
 * swap rows i and j in place
 */
void swaprows_inplace(matrix_t *m, int i, int j);

/**
 * swap columsn i and j in place
 */
void swapcols_inplace(matrix_t *m, int i, int j);



/**
 * take inverse of 3x3 matrix A and store it in matrix invA
 */
int inv3x3(matrix_t *A, matrix_t *invA);