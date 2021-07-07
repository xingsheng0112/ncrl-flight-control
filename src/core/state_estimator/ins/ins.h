#ifndef __INS_H__
#define __INS_H__

#include <stdbool.h>

void ins_init(void);
bool ins_check_sensor_status(void);
void ins_state_estimate(void);

/* raw position getters */
void ins_get_raw_position(float *pos);
float ins_get_raw_position_x(void);
float ins_get_raw_position_y(void);
float ins_get_raw_position_z(void);

/* raw velocity getters */
void ins_get_raw_velocity(float *vel);
float ins_get_raw_velocity_x(void);
float ins_get_raw_velocity_y(void);
float ins_get_raw_velocity_z(void);

/* ins fused position getters */
void ins_get_fused_position(float *pos);
float ins_get_fused_position_x(void);
float ins_get_fused_position_y(void);
float ins_get_fused_position_z(void);

/* ins fused velocity getters */
void ins_get_fused_velocity(float *vel);
float ins_get_fused_velocity_x(void);
float ins_get_fused_velocity_y(void);
float ins_get_fused_velocity_z(void);

/* ins ahrs attitude getters */
void ins_ahrs_get_attitude_euler_angles(float *roll, float *pitch, float *yaw);
void ins_ahrs_get_attitude_quaternion(float *q);
void ins_ahrs_get_rotation_matrix_b2i(float **R_b2i);
void ins_ahrs_get_rotation_matrix_i2b(float **R_i2b);

#endif
