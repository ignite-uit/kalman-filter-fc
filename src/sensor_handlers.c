#include <math.h>
#include "sensor_handlers.h"

static float sq(float x) {
    return x*x;
}

float barometer_altitude_variance(float pressure)
{
    //  R_g**2*T0**2/(g**2*M**2*p**2) * barometer_variance
    return R_g * R_g * T0 * T0 / (g * g * M * M * pressure * pressure) * barometer_variance;
}

float altitude(float pressure)
{
    return (R_g * T0) / (M * g) * (logf(P0 / pressure));
}

/**
 * Calculates the variance of the accelerometer measurement in the earth reference frame
   Calculation is based on the quaternion product

   sigma_ak must be allocated before use
 */
void ae_variance(quaternion_t *q, float *sigma_ak)
{
    float w = q->w, r1 = q->r1, r2 = q->r2, r3 = q->r3;

    float aex_variance = (4*sq(1/2.0f + w - sq(r2) - sq(r3)) + sq(r1) * sq(r2) + sq(r1) * sq(r3))*accelerometer_variance;
    float aey_variance = (4*sq(1/2.0f + w - sq(r1) - sq(r3)) + sq(r2) * sq(r3) + sq(r1) * sq(r2))*accelerometer_variance;
    float aez_variance = (4*sq(1/2.0f + w - sq(r1) - sq(r2)) + sq(r1) * sq(r3) + sq(r2) * sq(r3))*accelerometer_variance;

    // sigma_ak is statically allocated in sensor_handlers.h
    sigma_ak[0] = aex_variance;
    sigma_ak[1] = aey_variance;
    sigma_ak[2] = aez_variance;
    sigma_ak[3] = aex_variance;
    sigma_ak[4] = aey_variance;
    sigma_ak[5] = aez_variance;
}