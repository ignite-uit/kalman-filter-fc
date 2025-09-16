#include <stdio.h>
#include "sensor_handlers.h"
#include "kalman_config.h"
#include "math_util.h"

// static memory allocations for matrices
float F_data[dimState * dimState] = {0.0f}; // state model matrix
float B_data[dimState * numRowB] = {0.0f};  // control matrix
float H_data[numColH * numRowH] = {0.0f};   // observation matrix
float Q_data[dimState * dimState] = {0.0f}; // process noise matrix
float R_data[numRowR * numColR] = {0.0f};   // measurement noise matrix
float state_data[dimState] = {0.0f};
float sigma_ak[6] = {0.0f};

matrix_t F, B, H, Q, R;

static void updateR(float pressure);

void kalman_filter_init()
{
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
    F.numCol = dimState;
    F.numRow = dimState;
    F.data = F_data;
    for (int i = 0; i < dimState; i++)
        set_val(&F, i, i, 1.0);
    set_val(&F, 0, 3, Dt);
    set_val(&F, 1, 4, Dt);
    set_val(&F, 2, 5, Dt);

    // control matrix
    /*
          [1/2 * Dt**2, 0,           0          ],
          [0,           1/2 * Dt**2, 0          ],
     B =  [0,           0,           1/2 * Dt**2],
          [Dt,          0,           0          ],
          [0,           Dt,          0          ],
          [0,           0,           Dt         ]]
    */
    B.numCol = numColB;
    B.numRow = numRowB;
    B.data = B_data;
    for (int i = 0; i < numColB; i++)
        set_val(&B, i, i, 0.5f * Dt * Dt);
    set_val(&B, 3, 0, Dt);
    set_val(&B, 4, 1, Dt);
    set_val(&B, 5, 2, Dt);

    // observation matrix --- optimize if needed
    /*
             [1, 0, 0, 0, 0, 0],
       H =   [0, 1, 0, 0, 0, 0],
             [0, 0, 1, 0, 0, 0]]
    */
    H.numCol = numColH;
    H.numRow = numRowH;
    H.data = H_data;
    for (int i = 0; i < numColB; i++)
        set_val(&H, i, i, 1.0f);

    // Process noise matrix
    /*
                [1/4 * Dt**4, 0,           0,           1/2 * Dt**3, 0,           0          ],
                [0,           1/4 * Dt**4, 0,           0,           1/2 * Dt**3, 0          ],
         Q =    [0,           0,           1/4 * Dt**4, 0,           0,           1/2 * Dt**3],
                [1/2 * Dt**3, 0,           0,           Dt**2,       0,           0          ],
                [0,           1/2 * Dt**3, 0,           0,           Dt**2,       0          ],
                [0,           0,           1/2 * Dt**3, 0,           0,           Dt**2      ]]) * sigma_ak * Qgain
    */
    Q.numCol = dimState;
    Q.numRow = dimState;
    Q.data = Q_data;
    for (int i = 0; i < dimState; i++)
        set_val(&Q, i, i, 0.25f * Dt * Dt * Dt * Dt * accelerometer_variance * Qgain);
    set_val(&Q, 0, 3, 0.5f * Dt * Dt * Dt * accelerometer_variance * Qgain);
    set_val(&Q, 1, 4, 0.5f * Dt * Dt * Dt * accelerometer_variance * Qgain);
    set_val(&Q, 2, 5, 0.5f * Dt * Dt * Dt * accelerometer_variance * Qgain);

    set_val(&Q, 3, 0, 0.5f * Dt * Dt * Dt * accelerometer_variance * Qgain);
    set_val(&Q, 4, 1, 0.5f * Dt * Dt * Dt * accelerometer_variance * Qgain);
    set_val(&Q, 5, 2, 0.5f * Dt * Dt * Dt * accelerometer_variance * Qgain);

    /*
    [[GNSS_x_variance,  0,  0],
    [0, GNSS_y_variance,   0],
    [0,                0,  sigma_z]]) * Rgain    R.numCol = numColR;
    */

    R.numRow = numRowR;
    R.numCol = numColR;
    R.data = R_data;
    set_val(&R, 0, 0, GNSS_x_variance * Rgain);
    set_val(&R, 1, 1, GNSS_y_variance * Rgain);
    updateR(P0);
}

static void updateR(float pressure)
{
    // assumes the GNSS variance is constant.
    // Altitude variance changes with pressure
    float sigma_alt = barometer_altitude_variance(pressure);
    set_val(&R, 2, 2, sigma_alt * Rgain);
}



static void print_matrix(matrix_t *m)
{
    for (int i = 0; i < m->numRow; i++)
    {
        for (int j = 0; j < m->numCol; j++)
        {
            // printf("%f ", get_value(&result, i, j));
            printf("%f ", get_value(m, i, j));
        }
        printf("\n");
    }
}

static void clear_matrix_data(matrix_t *m) {
    for (int i = 0; i < m->numRow; i++)
    {
        for (int j = 0; j < m->numCol; j++)
        {
            set_val(m, i, j, 0.0f);
        }
    }

}

static void smoke_checks()
{
    // test that matmul works: R*H
    matrix_t result;
    result.numRow = H.numRow;
    result.numCol = R.numCol;
    float resultData[numRowH * numColR] = {0.0f};
    int errcode = 0;
    result.data = resultData;
    if (!matmul(&R, &H, &result, &errcode)) {
        printf("matrix multiplication failed");
    }

    printf("R matrix\n");
    print_matrix(&R);

    printf("H matrix\n");
    print_matrix(&H);

    printf("result matrix RxH\n");
    print_matrix(&result);
    printf("\n");

    clear_matrix_data(&result);
    errcode = 0;
    matmul(&H, &B, &result, &errcode);
    printf("B matrix\n");
    print_matrix(&B);
    printf("H x B\n");
    print_matrix(&result);

}

int main()
{
    kalman_filter_init();
    printf("hello kalman!\n");

    smoke_checks();

}
