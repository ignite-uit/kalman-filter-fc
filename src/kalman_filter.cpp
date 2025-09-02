

#define dimState 6
#define dimRowB 6
#define dimColB 3

#define dimColH 6
#define dimRowH 3

#define dimRowR 3
#define dimColR 3



#define sigma_ak 0.35*0.35 // double check
#define Qgain 1.0
#define Rgain 1.0

#define barometer_variance (0.5E2)*(0.5E2);

#define GNSS_x_variance 0.1
#define GNSS_y_variance 0.1
#define GNSS_xy_variance 0.1


#define g 9.81 // gravitational constant, unit: m/s^2
#define Dt 0.1 // timestep, unit: seconds
#define R_g 8.31432e3   // universal gas constant, unit: N m / (kmol K)
#define M 28.9644        // mean molecular weight of air at sea level, unit: kg/kmol

// Reference values, must be on the day before takeoff
#define P0 101325        // unit: Pa
#define T0 290           // unit: K (kelvin)

typedef struct matrix
{
    int dimCol;
    int dimRow;
    float *data;
} matrix_t;

typedef struct state_vector
{
    int dim;
    float *data;
} state_vector_t;

// static memory allocations for matrices
float F_data[dimState * dimState] = {0.0};  // state model matrix
float B_data[dimState * dimRowB] = {0.0};   // control matrix
float H_data[dimColH * dimRowH] = {0.0};    // observation matrix
float Q_data[dimState * dimState] = {0.0};  // process noise matrix
float R_data[dimRowR * dimColR] = {0.0};    // measurement noise matrix
float state_data[dimState] = {0.0};

matrix_t F, B, H, Q, R;

float barometer_altitude_variance(float pressure) {
    //  R_g**2*T0**2/(g**2*M**2*p**2) * barometer_variance
    return R_g*R_g*T0*T0/(g*g*M*M*pressure*pressure) * barometer_variance
}

void assign_value(matrix_t *m, int row, int col, float value) {
    int dimRow, dimCol, idx;
    dimRow = m->dimRow;
    dimCol = m->dimCol;

    idx = row*dimCol + col;

    m->data[idx] = value;
}


void kalman_filter_init() {
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
    F.dimCol = dimState;
    F.dimRow = dimState;
    F.data = F_data;
    for (int i=0; i<dimState; i++)
        assign_value(&F, i, i, 1.0);
    assign_value(&F, 0, 3, Dt);
    assign_value(&F, 1, 4, Dt);
    assign_value(&F, 2, 5, Dt);


    // control matrix
    /*
          [1/2 * Dt**2, 0,           0          ],
          [0,           1/2 * Dt**2, 0          ],
     B =  [0,           0,           1/2 * Dt**2],
          [Dt,          0,           0          ],
          [0,           Dt,          0          ],
          [0,           0,           Dt         ]]
    */
    B.dimCol = dimColB;
    B.dimRow = dimColB;
    B.data = B_data;
    for (int i=0; i<dimColB; i++)
        assign_value(&B, i, i, 0.5*Dt*Dt);
    assign_value(&B, 3, 0, Dt);
    assign_value(&B, 4, 1, Dt);
    assign_value(&B, 5, 2, Dt);


    // observation matrix --- optimize if needed
    /*
             [1, 0, 0, 0, 0, 0],
       H =   [0, 1, 0, 0, 0, 0],
             [0, 0, 1, 0, 0, 0]]  
    */
    H.dimCol = dimColH;
    H.dimRow = dimRowH;
    H.data = H_data;
    for (int i=0; i<dimColB; i++)
        assign_value(&B, i, i, 1.0);
    

    // Process noise matrix
    /*
                [1/4 * Dt**4, 0,           0,           1/2 * Dt**3, 0,           0          ],
                [0,           1/4 * Dt**4, 0,           0,           1/2 * Dt**3, 0          ],
         Q =    [0,           0,           1/4 * Dt**4, 0,           0,           1/2 * Dt**3],
                [1/2 * Dt**3, 0,           0,           Dt**2,       0,           0          ],
                [0,           1/2 * Dt**3, 0,           0,           Dt**2,       0          ],
                [0,           0,           1/2 * Dt**3, 0,           0,           Dt**2      ]]) * sigma_ak * Qgain 
    */
    Q.dimCol = dimState;
    Q.dimRow = dimState;
    Q.data = Q_data;
    for (int i=0; i < dimState; i++)
        assign_value(&Q, i, i, 0.25*Dt*Dt*Dt*Dt*sigma_ak*Qgain);
    assign_value(&Q, 0, 3, 0.5*Dt*Dt*Dt*sigma_ak*Qgain);
    assign_value(&Q, 1, 4, 0.5*Dt*Dt*Dt*sigma_ak*Qgain);
    assign_value(&Q, 2, 5, 0.5*Dt*Dt*Dt*sigma_ak*Qgain);

    assign_value(&Q, 3, 0, 0.5*Dt*Dt*Dt*sigma_ak*Qgain);
    assign_value(&Q, 4, 1, 0.5*Dt*Dt*Dt*sigma_ak*Qgain);
    assign_value(&Q, 5, 2, 0.5*Dt*Dt*Dt*sigma_ak*Qgain);

    R.dimCol = dimColR;
    R.dimRow = dimRowR;
    R.data = R_data;
    assign_value(&R, 0, 0, GNSS_x_variance);
    assign_value(&R, 1, 1, GNSS_y_variance);
    assign_value(&R, 0, 0, GNSS_x_variance);

}






int main() {
    kalman_filter_init();
}

