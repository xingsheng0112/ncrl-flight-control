#ifndef __MULTIROTOR_GEOMETRY_CTRL_H__
#define __MULTIROTOR_GEOMETRY_CTRL_H__

#include "imu.h"
#include "ahrs.h"
#include "debug_link.h"
#include "matrix.h"

#define ICL_N 20

typedef struct ICL_sigma{
	float y;
	MAT_ALLOC(y_omega,3,1);	
	MAT_ALLOC(M_hat,3,1);	
} ICL_sigma_t;

void geometry_ctrl_init(void);
void multirotor_geometry_control(radio_t *rc, float *desired_heading);

void send_geometry_moment_ctrl_debug(debug_msg_t *payload);
void send_geometry_tracking_ctrl_debug(debug_msg_t *payload);
void send_uav_dynamics_debug(debug_msg_t *payload);
void send_controller_estimation_adaptive_debug(debug_msg_t *payload);

#endif
