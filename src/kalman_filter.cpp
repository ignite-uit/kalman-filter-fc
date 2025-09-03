#include <stdio.h>
#include "sensor_handlers.h"
#include "kalman_config.h"
#include "math_util.h"

// static memory allocations for matrices
float F_data[dimState * dimState] = {0.0f}; // state model matrix
float B_data[dimState * dimRowB] = {0.0f};  // control matrix
float H_data[dimColH * dimRowH] = {0.0f};   // observation matrix
float Q_data[dimState * dimState] = {0.0f}; // process noise matrix
float R_data[dimRowR * dimColR] = {0.0f};   // measurement noise matrix
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
    F.dimCol = dimState;
    F.dimRow = dimState;
    F.data = F_data;
    for (int i = 0; i < dimState; i++)
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
    B.dimRow = dimRowB;
    B.data = B_data;
    for (int i = 0; i < dimColB; i++)
        assign_value(&B, i, i, 0.5f * Dt * Dt);
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
    for (int i = 0; i < dimColB; i++)
        assign_value(&H, i, i, 1.0f);

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
    for (int i = 0; i < dimState; i++)
        assign_value(&Q, i, i, 0.25f * Dt * Dt * Dt * Dt * accelerometer_variance * Qgain);
    assign_value(&Q, 0, 3, 0.5f * Dt * Dt * Dt * accelerometer_variance * Qgain);
    assign_value(&Q, 1, 4, 0.5f * Dt * Dt * Dt * accelerometer_variance * Qgain);
    assign_value(&Q, 2, 5, 0.5f * Dt * Dt * Dt * accelerometer_variance * Qgain);

    assign_value(&Q, 3, 0, 0.5f * Dt * Dt * Dt * accelerometer_variance * Qgain);
    assign_value(&Q, 4, 1, 0.5f * Dt * Dt * Dt * accelerometer_variance * Qgain);
    assign_value(&Q, 5, 2, 0.5f * Dt * Dt * Dt * accelerometer_variance * Qgain);

    /*
    [[GNSS_x_variance,  0,  0],
    [0, GNSS_y_variance,   0],
    [0,                0,  sigma_z]]) * Rgain    R.dimCol = dimColR;
    */

    R.dimRow = dimRowR;
    R.dimCol = dimColR;
    R.data = R_data;
    assign_value(&R, 0, 0, GNSS_x_variance * Rgain);
    assign_value(&R, 1, 1, GNSS_y_variance * Rgain);
    updateR(P0);
}

static void updateR(float pressure)
{
    // assumes the GNSS variance is constant.
    // Altitude variance changes with pressure
    float sigma_alt = barometer_altitude_variance(pressure);
    assign_value(&R, 2, 2, sigma_alt * Rgain);
}

static void print_matrix(matrix_t *m)
{
    for (int i = 0; i < m->dimRow; i++)
    {
        for (int j = 0; j < m->dimCol; j++)
        {
            // printf("%f ", get_value(&result, i, j));
            printf("%f ", get_value(m, i, j));
        }
        printf("\n");
    }
}

static void clear_matrix_data(matrix_t *m) {
    for (int i = 0; i < m->dimRow; i++)
    {
        for (int j = 0; j < m->dimCol; j++)
        {
            assign_value(m, i, j, 0.0f);
        }
    }

}

static void smoke_checks()
{
    // test that matmul works: R*H
    matrix_t result;
    result.dimRow = H.dimRow;
    result.dimCol = R.dimCol;
    float resultData[dimRowH * dimColR] = {0.0f};
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
