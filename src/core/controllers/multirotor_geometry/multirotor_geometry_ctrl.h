#ifndef __MULTIROTOR_GEOMETRY_CTRL_H__
#define __MULTIROTOR_GEOMETRY_CTRL_H__

#include "imu.h"
#include "ahrs.h"
#include "debug_link.h"
#include "matrix.h"

#define ICL_N 10
//25 for time step

typedef struct {
	MAT_ALLOC(y_cl,4,6);	
	MAT_ALLOC(F,4,1);	//F_bar
} ICL_sigma;

void geometry_ctrl_init(void);
void multirotor_geometry_control(radio_t *rc);

void send_geometry_moment_ctrl_debug(debug_msg_t *payload);
void send_geometry_tracking_ctrl_debug(debug_msg_t *payload);
void send_uav_dynamics_debug(debug_msg_t *payload);
void send_controller_estimation_adaptive_debug(debug_msg_t *payload);

#endif
