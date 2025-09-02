

#define dimState 6
#define dimColH 3
#define dimRowB 3

#define dimRowR 3
#define dimColR 3

#define dimRowH 6


typedef struct matrix {
    int dimCol;
    int dimRow;
    float *data;
} matrix_t;

typedef struct state_vector {
    int dim;
    float *data;
} state_vector_t;


float F_data[dimState*dimState];
float B_data[dimState*dimRowB];
float Q_data[dimState*dimState];
float R_data[dimRowR*dimColR];

matrix_t F, B, H, Q, R;

float Dt = 0.1; // timestep


// static allocations of kalman filter matrices
/*
state model matrix
Maps the previous state vector to the next during the prediction step
    [[1, 0, 0, Dt, 0,  0 ],
     [0, 1, 0, 0,  Dt, 0 ],
     [0, 0, 1, 0,  0,  Dt],
     [0, 0, 0, 1,  0,  0 ],
     [0, 0, 0, 0,  1,  0 ],
     [0, 0, 0, 0,  0,  1 ]]
*/




