
#define accelerometer_variance 0.35f * 0.35f

#define barometer_variance (0.5E2f) * (0.5E2f)

#define GNSS_x_variance 0.1f
#define GNSS_y_variance 0.1f
#define GNSS_xy_variance 0.1f

#define g 9.81f        // gravitational constant, unit: m/s^2
#define Dt 0.1f        // timestep, unit: seconds
#define R_g 8.31432e3f // universal gas constant, unit: N m / (kmol K)
#define M 28.9644f     // mean molecular weight of air at sea level, unit: kg/kmol

// Reference values, must be on the day before takeoff
#define P0 101325.0f // unit: Pa
#define T0 290.0f    // unit: K (kelvin)


typedef struct quaternion {
    float w;
    float r1;
    float r2;
    float r3;
} quaternion_t;

float barometer_altitude_variance(float pressure);