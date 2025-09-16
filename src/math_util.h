
#define MATMUL_DIMENSION_MISMATCH_ERROR 1
#define MATADD_DIMENSION_MISMATCH_ERROR 2
#define MAT_INV_SINGULAR_MATRIX_ERROR 3

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


/**
 * Generates a variable prefix_name
 * Example usage: int GENERATE_VAR(my, value) = 5; // similar to int my_value = 5;
 */
#define GENERATE_VAR(prefix, name) prefix ## _ ## name

/**
 * Allocate a vector_t vectorname with dimension dim=dim and its data float[dim] statically
 */
#define stackVectorAllocate(name, dim) \
    float GENERATE_VAR(name, data)[dim]; \
    vector_t name; \
    (name).data = GENERATE_VAR(name, data);

#define stackMatrixAllocate(name, numRow, numCol) \
    float GENERATE_VAR(name, data)[numRow * numCol]; \
    matrix_t name; \
    (name).data = GENERATE_VAR(name, data);

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

void mult_mat_scal(matrix_t *mat, float scalar);

int matsub(matrix_t *a, matrix_t *b, matrix_t *result, int *errorcode);
int matadd(matrix_t *a, matrix_t *b, matrix_t *result, int *errorcode);
int vecadd(vector_t *a, vector_t *b, vector_t *output, int *errorcode);
int vecsub(vector_t *a, vector_t *b, vector_t *output, int *errorcode);
int matvecmul(matrix_t *A, vector_t *input, vector_t *output, int *errorcode);



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
int inv3x3(matrix_t *A, matrix_t *invA, int *errorcode);