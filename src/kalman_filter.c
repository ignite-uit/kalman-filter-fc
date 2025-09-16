#include <stdio.h>
#include "sensor_handlers.h"
#include "kalman_config.h"
#include "math_util.h"

// static memory allocations for matrices
float Id_data[dimState * dimState] = {0.0f}; // Identity matrix
float F_data[dimState * dimState] = {0.0f};  // state model matrix
float Ft_data[dimState * dimState] = {0.0f}; // state model matrix
float B_data[dimState * numRowB] = {0.0f};   // control matrix
float H_data[numColH * numRowH] = {0.0f};    // observation matrix
float Ht_data[numColH * numRowH] = {0.0f};   // observation matrix transpose
float Q_data[dimState * dimState] = {0.0f};  // process noise matrix
float R_data[numRowR * numColR] = {0.0f};    // measurement noise matrix
float P_data[dimState * dimState] = {1.0f};  // prediction covariance matrix
float state_data[dimState] = {0.0f};
float sigma_ak[6] = {0.0f};

matrix_t Id, F, Ft, B, H, Ht, Q, R, P;

// statically allocate space for the residuals
float yk_data[dimState] = {0.0f};
vector_t yk;

// statically allocate space for the state vector
float xkk_data[dimState] = {0.0f};
vector_t xkk;

static void updateR(float pressure)
{
    // assumes the GNSS variance is constant.
    // Altitude variance changes with pressure
    // sigma_zk = [sigma_x, sigma_y, sigma_alt]
    float sigma_alt = barometer_altitude_variance(pressure);
    set_val(&R, 0, 0, GNSS_x_variance);
    set_val(&R, 1, 1, GNSS_y_variance);
    set_val(&R, 2, 2, sigma_alt * Rgain);
}

void kalman_filter_init()
{

    yk.dim = dimState;
    yk.data = yk_data;
    xkk.dim = dimState;
    xkk.data = xkk_data;

    Id.numCol = dimState;
    Id.numRow = dimState;
    Id.data = Id_data;
    for (int i = 0; i < dimState; i++)
    {
        set_val(&Id, i, i, 1.0f);
    }

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
    {
        set_val(&F, i, i, 1.0);
        set_val(&F, 0, 3, Dt);
        set_val(&F, 1, 4, Dt);
        set_val(&F, 2, 5, Dt);
    }

    Ft.numCol = dimState;
    Ft.numRow = dimState;
    Ft.data = F_data;

    for (int i = 0; i < dimState; i++)
    {
        set_val(&Ft, i, i, 1.0);
        set_val(&Ft, 3, 0, Dt);
        set_val(&Ft, 4, 1, Dt);
        set_val(&Ft, 5, 2, Dt);
    }

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
    {
        set_val(&B, i, i, 0.5f * Dt * Dt);
        set_val(&B, 3, 0, Dt);
        set_val(&B, 4, 1, Dt);
        set_val(&B, 5, 2, Dt);
    }

    // observation matrix --- optimize if needed
    /*
             [1, 0, 0, 0, 0, 0],
       H =   [0, 1, 0, 0, 0, 0],
             [0, 0, 1, 0, 0, 0]]
    */
    H.numCol = numColH;
    H.numRow = numRowH;
    H.data = H_data;
    for (int i = 0; i < numColH; i++)
        set_val(&H, i, i, 1.0f);

    // transpose observation matrix --- optimize if needed
    /*
             [1, 0, 0],
       H =   [0, 1, 0],
             [0, 0, 1]
             [0, 0, 0]
             [0, 0, 0]
             [0, 0, 0]
    */
    Ht.numCol = numColH;
    Ht.numRow = numRowH;
    Ht.data = H_data;
    for (int i = 0; i < numColH; i++)
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
    {
        set_val(&Q, i, i, 0.25f * Dt * Dt * Dt * Dt * accelerometer_variance * Qgain);
        set_val(&Q, 0, 3, 0.5f * Dt * Dt * Dt * accelerometer_variance * Qgain);
        set_val(&Q, 1, 4, 0.5f * Dt * Dt * Dt * accelerometer_variance * Qgain);
        set_val(&Q, 2, 5, 0.5f * Dt * Dt * Dt * accelerometer_variance * Qgain);

        set_val(&Q, 3, 0, 0.5f * Dt * Dt * Dt * accelerometer_variance * Qgain);
        set_val(&Q, 4, 1, 0.5f * Dt * Dt * Dt * accelerometer_variance * Qgain);
        set_val(&Q, 5, 2, 0.5f * Dt * Dt * Dt * accelerometer_variance * Qgain);
    }

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

    P.numRow = dimState;
    P.numCol = dimState;
}

int getQgain(float *qgain)
{
    *qgain = Qgain;
    return 0;
}

int predict(vector_t *predVec, matrix_t *predCov, vector_t *ak, int *errorcode)
{
    stackVectorAllocate(Fx_k, dimState);
    stackVectorAllocate(Ba_k, dimState);
    stackMatrixAllocate(FP, dimState, dimState);
    stackMatrixAllocate(FPFt, dimState, dimState);

    // update prediciton vector predVec
    ////////////////////////////////////////////////
    // F * state vector (our global xkk)
    if (!matvecmul(&F, &xkk, &Fx_k, errorcode))
        goto cleanup;

    // Ba_k = B*ak
    if (!matvecmul(&B, ak, &Ba_k, errorcode))
        goto cleanup;

    // predVec = F * statevec - Ba_k
    if (!vecsub(&Fx_k, &Ba_k, predVec, errorcode))
        goto cleanup;
    ////////////////////////////////////////////////

    // update prediction matrix predCov using our current
    // covariance matrix P using F and ading Q
    // predCov = F * P_current * F.T + Q
    if (!matmul(&F, &P, &FP, errorcode))
        goto cleanup;

    if (!matmul(&FP, &Ft, &FPFt, errorcode))
        goto cleanup;

    if (!matadd(&FPFt, &Q, predCov, errorcode))
        goto cleanup;

    return 1;

cleanup:
    fprintf(stderr, "kalman: predict failed, errorcode: %d\n", *errorcode);
    return 0;
}

int update(vector_t *predVec, matrix_t *pred_cov_mat, vector_t *zk, float pressure, int *errorcode)
{
    stackVectorAllocate(Hx, numRowH);
    stackVectorAllocate(Ky_k, dimState);

    // residual covariances
    stackMatrixAllocate(Sk, numRowR, numColR);
    stackMatrixAllocate(invSk, numRowR, numColR);

    // Kalman gain
    stackMatrixAllocate(Kk, dimState, numRowH);
    stackMatrixAllocate(Kk_times_H, dimState, numColH); // dimState x dimState
    stackMatrixAllocate(Id_minus_Kk_times_H, dimState, numColH);

    stackMatrixAllocate(Pkkm1_X_Ht, pred_cov_mat->numCol, numColH);
    stackMatrixAllocate(H_times_Pkkm1_times_Ht, pred_cov_mat->numCol, numColH);
    stackMatrixAllocate(Ht_X_invSK, numColH, numColR);

    updateR(pressure);

    // yk = zk - H * predVec
    if (!matvecmul(&H, predVec, &Hx, errorcode))
        goto errorcleanup;

    // store yk = zk - Hx
    if (!vecsub(zk, &Hx, &yk, errorcode))
        goto errorcleanup;

    // calculate residual covariance
    /////////////////////////////////////////////////////
    // Sk = H * Pkkm1 * H.t + R
    if (!matmul(pred_cov_mat, &Ht, &Pkkm1_X_Ht, errorcode))
        goto errorcleanup;
    // left multiply by H
    if (!matmul(&H, &Pkkm1_X_Ht, &H_times_Pkkm1_times_Ht, errorcode))
        goto errorcleanup;
    // add measurement covariance matrix R
    if (!matadd(&H_times_Pkkm1_times_Ht, &R, &Sk, errorcode))
        goto errorcleanup;
    /////////////////////////////////////////////////////

    // Calculate Kalman gain
    ////////////////////////////////////////////////////
    // calculate inverse
    if (!inv3x3(&Sk, &invSk, errorcode))
        goto errorcleanup;

    // H.t * inv(Sk)
    if (!matmul(&Ht, &invSk, &Ht_X_invSK, errorcode))
        goto errorcleanup;

    if (!matmul(pred_cov_mat, &Ht_X_invSK, &Kk, errorcode))
        goto errorcleanup;
    ////////////////////////////////////////////////////

    // calculate new estimate
    /////////////////////////////////////////////////////
    // calculate Ky_k = Kk * yk
    if (!matvecmul(&Kk, &yk, &Ky_k, errorcode))
        goto errorcleanup;

    // xkk = predVec + Ky_k
    if (!vecadd(predVec, &Ky_k, &xkk, errorcode))
        goto errorcleanup;
    /////////////////////////////////////////////////////

    // Calculate new prediction covariance matrix
    // Pkk = (Id - (Kk * H))* pred_cov_mat
    /////////////////////////////////////////////////////
    // Kk * H
    if (!matmul(&Kk, &H, &Kk_times_H, errorcode))
        goto errorcleanup;

    if (!matsub(&Id, &Kk_times_H, &Id_minus_Kk_times_H, errorcode))
        goto errorcleanup;

    if (!matmul(&Id_minus_Kk_times_H, pred_cov_mat, &P, errorcode))
        goto errorcleanup;
    /////////////////////////////////////////////////////

    // update residuals
    // yk = zk - H*xkk
    if (!matvecmul(&H, &xkk, &Hx, errorcode))
        goto errorcleanup;

    return 1;

errorcleanup:
    fprintf(stderr, "kalman update function failed,"
                    "errorcode %d, \n",
            *errorcode);
    return 0;
}

int KF_one_iteration(vector_t *ak, vector_t *zk, float pressure, int *errorcode) {
    // ak -- accelerometer data in meters per second and in earth frame of reference
    // zk -- k'th GNSS and barometer measurement in meters
    //  int predict(vector_t *predVec, matrix_t *predCov, vector_t *ak, int *errorcode)
    stackVectorAllocate(pred_vec, dimState);
    stackMatrixAllocate(pred_cov_mat, dimState, dimState);
    if (!predict(&pred_vec, &pred_cov_mat, ak, errorcode))
        return 1;
    if (!update(&pred_vec, &pred_cov_mat, zk, pressure,  errorcode))
        return 1;
    return 0;
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
    if (!matmul(&R, &H, &result, &errcode))
    {
        printf("matrix multiplication failed");
    }

    printf("R matrix\n");
    pprint_matrix(&R);

    printf("H matrix\n");
    pprint_matrix(&H);

    printf("result matrix RxH\n");
    pprint_matrix(&result);
    printf("\n");

    clear_matrix(&result);
    errcode = 0;
    matmul(&H, &B, &result, &errcode);
    printf("B matrix\n");
    pprint_matrix(&B);
    printf("H x B\n");
    pprint_matrix(&result);
}

int main()
{
    kalman_filter_init();
    printf("hello kalman!\n");

    smoke_checks();
}
