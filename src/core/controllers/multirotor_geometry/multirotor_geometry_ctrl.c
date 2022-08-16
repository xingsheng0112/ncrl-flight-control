#include "proj_config.h"
#include "arm_math.h"
#include "FreeRTOS.h"
#include "task.h"
#include "semphr.h"
#include "gpio.h"
#include "sbus_radio.h"
#include "ahrs.h"
#include "matrix.h"
#include "motor_thrust_fitting.h"
#include "motor.h"
#include "bound.h"
#include "se3_math.h"
#include "lpf.h"
#include "imu.h"
#include "ahrs.h"
#include "autopilot.h"
#include "debug_link.h"
#include "multirotor_geometry_param.h"
#include "multirotor_geometry_ctrl.h"
#include "position_state.h"
#include "multirotor_rc.h"
#include "barometer.h"
#include "compass.h"
#include "sys_param.h"
#include "led.h"
#include "attitude_state.h"
#include "waypoint_following.h"
#include "fence.h"


#define dt 0.0025 //[s]
#define MOTOR_TO_CG_LENGTH 16.25f //[cm]
#define MOTOR_TO_CG_LENGTH_M (MOTOR_TO_CG_LENGTH * 0.01) //[m]
#define COEFFICIENT_YAW 0.117f

MAT_ALLOC(J, 3, 3);
MAT_ALLOC(R, 3, 3);
MAT_ALLOC(Rd, 3, 3);
MAT_ALLOC(Rt, 3, 3);
MAT_ALLOC(Rtd, 3, 3);
MAT_ALLOC(RtdR, 3, 3);
MAT_ALLOC(RtRd, 3, 3);
MAT_ALLOC(RtRdWd, 3, 3);
MAT_ALLOC(Re3, 3, 1);
MAT_ALLOC(W, 3, 1);
MAT_ALLOC(W_dot, 3, 1);
MAT_ALLOC(Wd, 3, 1);
MAT_ALLOC(W_hat, 3, 3);
MAT_ALLOC(Wd_dot, 3, 1);
MAT_ALLOC(JW, 3, 1);
MAT_ALLOC(WJW, 3, 1);
MAT_ALLOC(JWdot, 3, 1);
MAT_ALLOC(M, 3, 1);
MAT_ALLOC(eR_mat, 3, 3);
MAT_ALLOC(eR, 3, 1);
MAT_ALLOC(eW, 3, 1);
MAT_ALLOC(WRt, 3, 3);
MAT_ALLOC(WRtRd, 3, 3);
MAT_ALLOC(WRtRdWd, 3, 1);
MAT_ALLOC(RtRdWddot, 3, 1);
MAT_ALLOC(WRtRdWd_RtRdWddot, 3, 1);
MAT_ALLOC(J_WRtRdWd_RtRdWddot, 3, 1);
MAT_ALLOC(inertia_effect, 3, 1);
MAT_ALLOC(kxex_kvev_mge3_mxd_dot_dot, 3, 1);
MAT_ALLOC(b1d, 3, 1);
MAT_ALLOC(b2d, 3, 1);
MAT_ALLOC(b3d, 3, 1);

float pos_error[3];
float vel_error[3];
float tracking_error_integral[3];

float krx, kry, krz;
float kwx, kwy, kwz;
float kpx, kpy, kpz;
float kvx, kvy, kvz;
float yaw_rate_ctrl_gain;
float k_tracking_i_gain[3];

float uav_mass;

//M = (J * W_dot) + (W X JW)
float uav_dynamics_m[3] = {0.0f};
//M_rot = (J * W_dot)
float uav_dynamics_m_rot_frame[3] = {0.0f};

float coeff_cmd_to_thrust[6] = {0.0f};
float coeff_thrust_to_cmd[6] = {0.0f};
float motor_thrust_max = 0.0f;

bool height_ctrl_only = false;

//adaptive & ICL parameters

//adaptive
MAT_ALLOC(theta, 8, 1);
MAT_ALLOC(theta_hat_dot, 8, 1);
MAT_ALLOC(Y, 3, 8);
MAT_ALLOC(Yt, 8, 3);
MAT_ALLOC(eW_c2eR, 3, 1);
MAT_ALLOC(Yt_eW_c2eR,8,1);
MAT_ALLOC(Ytheta,3,1);
//float adaptive_gamma[8] = {0.000005,0.000005,0.000028,0.000028,0.000005,0.0000028,0.0000028,0.0000028};
float adaptive_gamma[8] = {0.000002,0.000002,0.000002,0.000002,0.000002,0.000002,0.000002,0.000002};

float c2 = 10;

//ICL
MAT_ALLOC(adaptive_theta_hat_dot, 8, 1);
MAT_ALLOC(ICL_theta_hat_dot, 8, 1);
MAT_ALLOC(W_last,3,1);
MAT_ALLOC(W_last2,3,1);
MAT_ALLOC(W_last3,3,1);
MAT_ALLOC(W_last4,3,1);
MAT_ALLOC(W_last5,3,1);
MAT_ALLOC(W_dot,3,1);
MAT_ALLOC(Y1_Y_W, 3, 8);
MAT_ALLOC(Y1_Y_Wdt,3,8);
MAT_ALLOC(W_dot_Matrix, 3, 8);
MAT_ALLOC(M_last,3,1);
MAT_ALLOC(M_last2,3,1);
MAT_ALLOC(M_last3,3,1);
MAT_ALLOC(M_last4,3,1);
MAT_ALLOC(M_last5,3,1);
MAT_ALLOC(Mb,3,1);
MAT_ALLOC(y_cl,3,8);
//sigma
MAT_ALLOC(ICL_control_term, 8, 1);	//x
MAT_ALLOC(y_cltheta,3,1);			//y_cl*theta
MAT_ALLOC(Mb_y_cltheta,3,1);		//Mb-y_cl*theta
MAT_ALLOC(y_clT,8,3);				//y_cl.transpose
MAT_ALLOC(y_clT_Mb_y_cltheta,8,1);	//y_clT*Mb-y_cl*theta

float output_force_last = 0;
float output_force_last1 = 0;
float output_force_last2 = 0;
float output_force_last3 = 0;
float output_force_last4 = 0;
float output_force_last5 = 0;
//float k_icl[8] = {0.1,0.1,1200,1200,2200, 4000, 4000, 4000};
float k_icl[8] = {20000,20000,20000,20000,20000, 20000,20000,20000};
float adaptive_gamma_k_icl[8];


int ICL_sigma_index = 0;
ICL_sigma sigma_array[ICL_N];


void geometry_ctrl_init(void)
{
	init_multirotor_geometry_param_list();

	autopilot_init();

	float geo_fence_origin[3] = {0.0f, 0.0f, 0.0f};
	autopilot_set_enu_rectangular_fence(geo_fence_origin, 2.5f, 1.3f, 3.0f);

	MAT_INIT(J, 3, 3);
	MAT_INIT(R, 3, 3);
	MAT_INIT(Rd, 3, 3);
	MAT_INIT(Rt, 3, 3);
	MAT_INIT(Rtd, 3, 3);
	MAT_INIT(RtdR, 3, 3);
	MAT_INIT(RtRd, 3, 3);
	MAT_INIT(RtRdWd, 3, 3);
	MAT_INIT(Re3, 3, 1);
	MAT_INIT(W, 3, 1);
	MAT_INIT(W_dot, 3, 1);
	MAT_INIT(Wd, 3, 1);
	MAT_INIT(W_hat, 3, 3);
	MAT_INIT(Wd_dot, 3, 1);
	MAT_INIT(JW, 3, 1);
	MAT_INIT(WJW, 3, 1);
	MAT_INIT(JWdot, 3, 1);
	MAT_INIT(M, 3, 1);
	MAT_INIT(eR_mat, 3, 3);
	MAT_INIT(eR, 3, 1);
	MAT_INIT(eW, 3, 1);
	MAT_INIT(WRt, 3, 3);
	MAT_INIT(WRtRd, 3, 3);
	MAT_INIT(WRtRdWd, 3, 1);
	MAT_INIT(RtRdWddot, 3, 1);
	MAT_INIT(WRtRdWd_RtRdWddot, 3, 1);
	MAT_INIT(J_WRtRdWd_RtRdWddot, 3, 1);
	MAT_INIT(inertia_effect, 3, 1);
	MAT_INIT(kxex_kvev_mge3_mxd_dot_dot, 3, 1);
	MAT_INIT(b1d, 3, 1);
	MAT_INIT(b2d, 3, 1);
	MAT_INIT(b3d, 3, 1);

//adaptive & ICL parameters
//adaptive
	MAT_INIT(theta, 8, 1);
	mat_data(theta)[0] = 0.0f;
	mat_data(theta)[1] = 0.0f;
	mat_data(theta)[2] = 0.006f;
	mat_data(theta)[3] = 0.006f;
	mat_data(theta)[4] = 0.011f;
	mat_data(theta)[5] = 0.0f;
	mat_data(theta)[6] = 0.0f;
	mat_data(theta)[7] = 0.0f;
	MAT_INIT(theta_hat_dot, 8, 1);
	MAT_INIT(Y, 3, 8);
	MAT_INIT(Yt, 8, 3);
	MAT_INIT(eW_c2eR, 3, 1);
	MAT_INIT(Yt_eW_c2eR,8,1);
	MAT_INIT(Ytheta,3,1);
//ICL
	MAT_INIT(adaptive_theta_hat_dot, 8, 1);
	MAT_INIT(ICL_theta_hat_dot, 8, 1);
	MAT_INIT(W_last,3,1);
	MAT_INIT(W_last2,3,1);
	MAT_INIT(W_last3,3,1);
	MAT_INIT(W_last4,3,1);
	MAT_INIT(W_last5,3,1);
	mat_data(W_last)[0] = 0;
	mat_data(W_last)[1] = 0;
	mat_data(W_last)[2] = 0;
	mat_data(W_last2)[0] = 0;
	mat_data(W_last2)[1] = 0;
	mat_data(W_last2)[2] = 0;
	mat_data(W_last3)[0] = 0;
	mat_data(W_last3)[1] = 0;
	mat_data(W_last3)[2] = 0;
	mat_data(W_last4)[0] = 0;
	mat_data(W_last4)[1] = 0;
	mat_data(W_last4)[2] = 0;
	mat_data(W_last5)[0] = 0;
	mat_data(W_last5)[1] = 0;
	mat_data(W_last5)[2] = 0;
	MAT_INIT(Y1_Y_W, 3, 8);
	MAT_INIT(Y1_Y_Wdt,3,8);
	MAT_INIT(W_dot,3,1);
	MAT_INIT(W_dot_Matrix, 3, 8);
	MAT_INIT(M_last,3,1);
	MAT_INIT(M_last2,3,1);
	MAT_INIT(M_last3,3,1);
	MAT_INIT(M_last4,3,1);
	MAT_INIT(M_last5,3,1);
	mat_data(M_last)[0] = 0;
	mat_data(M_last)[1] = 0;
	mat_data(M_last)[2] = 0;
	mat_data(M_last2)[0] = 0;
	mat_data(M_last2)[1] = 0;
	mat_data(M_last2)[2] = 0;
	mat_data(M_last3)[0] = 0;
	mat_data(M_last3)[1] = 0;
	mat_data(M_last3)[2] = 0;
	mat_data(M_last4)[0] = 0;
	mat_data(M_last4)[1] = 0;
	mat_data(M_last4)[2] = 0;
	mat_data(M_last5)[0] = 0;
	mat_data(M_last5)[1] = 0;
	mat_data(M_last5)[2] = 0;
	MAT_INIT(Mb,3,1);
	MAT_INIT(y_cl,3,8);
	MAT_INIT(ICL_control_term, 8, 1);	//x
	MAT_INIT(y_cltheta,3,1);			//y_cl*theta
	MAT_INIT(Mb_y_cltheta,3,1);			//Mb-y_cl*theta
	MAT_INIT(y_clT,8,3);				//y_cl.transpose
	MAT_INIT(y_clT_Mb_y_cltheta,8,1);	//y_clT*Mb-y_cl*theta

	for( int i = 0 ; i<ICL_N ; i++ ){
			MAT_INIT(sigma_array[i].y_cl,3,8);
			MAT_INIT(sigma_array[i].Mb,3,1);
			
			for( int j = 0 ; j<3*8 ; j++){
				sigma_array[i].mat_data(y_cl)[j] = 0;
			}	
			sigma_array[i].mat_data(Mb)[0] = 0;
			sigma_array[i].mat_data(Mb)[1] = 0;
			sigma_array[i].mat_data(Mb)[2] = 0;
	}
	for( int i = 0 ; i<8 ; i++ ){
		adaptive_gamma_k_icl[i] = adaptive_gamma[i]*k_icl[i];
	}
	/* modify local variables when user change them via ground station */
	set_sys_param_update_var_addr(MR_GEO_GAIN_ROLL_P, &krx);
	set_sys_param_update_var_addr(MR_GEO_GAIN_ROLL_D, &kwx);
	set_sys_param_update_var_addr(MR_GEO_GAIN_PITCH_P, &kry);
	set_sys_param_update_var_addr(MR_GEO_GAIN_PITCH_D, &kwy);
	set_sys_param_update_var_addr(MR_GEO_GAIN_YAW_P, &krz);
	set_sys_param_update_var_addr(MR_GEO_GAIN_YAW_D, &kwz);
	set_sys_param_update_var_addr(MR_GEO_GAIN_RATE_YAW, &yaw_rate_ctrl_gain);
	set_sys_param_update_var_addr(MR_GEO_GAIN_POS_X, &kpx);
	set_sys_param_update_var_addr(MR_GEO_GAIN_VEL_X, &kvx);
	set_sys_param_update_var_addr(MR_GEO_GAIN_POS_Y, &kpy);
	set_sys_param_update_var_addr(MR_GEO_GAIN_VEL_Y, &kvy);
	set_sys_param_update_var_addr(MR_GEO_GAIN_POS_Z, &kpz);
	set_sys_param_update_var_addr(MR_GEO_GAIN_VEL_Z, &kvz);
	set_sys_param_update_var_addr(MR_GEO_GAIN_POS_X_I, &k_tracking_i_gain[0]);
	set_sys_param_update_var_addr(MR_GEO_GAIN_POS_Y_I, &k_tracking_i_gain[1]);
	set_sys_param_update_var_addr(MR_GEO_GAIN_POS_Z_I, &k_tracking_i_gain[2]);
	set_sys_param_update_var_addr(MR_GEO_UAV_MASS, &uav_mass);
	set_sys_param_update_var_addr(MR_GEO_INERTIA_JXX, &mat_data(J)[0*3 + 0]);
	set_sys_param_update_var_addr(MR_GEO_INERTIA_JYY, &mat_data(J)[1*3 + 1]);
	set_sys_param_update_var_addr(MR_GEO_INERTIA_JZZ, &mat_data(J)[2*3 + 2]);
	set_sys_param_update_var_addr(PWM_TO_THRUST_C1, &coeff_cmd_to_thrust[0]);
	set_sys_param_update_var_addr(PWM_TO_THRUST_C2, &coeff_cmd_to_thrust[1]);
	set_sys_param_update_var_addr(PWM_TO_THRUST_C3, &coeff_cmd_to_thrust[2]);
	set_sys_param_update_var_addr(PWM_TO_THRUST_C4, &coeff_cmd_to_thrust[3]);
	set_sys_param_update_var_addr(PWM_TO_THRUST_C5, &coeff_cmd_to_thrust[4]);
	set_sys_param_update_var_addr(PWM_TO_THRUST_C6, &coeff_cmd_to_thrust[5]);
	set_sys_param_update_var_addr(THRUST_TO_PWM_C1, &coeff_thrust_to_cmd[0]);
	set_sys_param_update_var_addr(THRUST_TO_PWM_C2, &coeff_thrust_to_cmd[1]);
	set_sys_param_update_var_addr(THRUST_TO_PWM_C3, &coeff_thrust_to_cmd[2]);
	set_sys_param_update_var_addr(THRUST_TO_PWM_C4, &coeff_thrust_to_cmd[3]);
	set_sys_param_update_var_addr(THRUST_TO_PWM_C5, &coeff_thrust_to_cmd[4]);
	set_sys_param_update_var_addr(THRUST_TO_PWM_C6, &coeff_thrust_to_cmd[5]);
	set_sys_param_update_var_addr(THRUST_MAX, &motor_thrust_max);

	/* load local variables previously stored in internal flash */
	get_sys_param_float(MR_GEO_GAIN_ROLL_P, &krx);
	get_sys_param_float(MR_GEO_GAIN_ROLL_D, &kwx);
	get_sys_param_float(MR_GEO_GAIN_PITCH_P, &kry);
	get_sys_param_float(MR_GEO_GAIN_PITCH_D, &kwy);
	get_sys_param_float(MR_GEO_GAIN_YAW_P, &krz);
	get_sys_param_float(MR_GEO_GAIN_YAW_D, &kwz);
	get_sys_param_float(MR_GEO_GAIN_RATE_YAW, &yaw_rate_ctrl_gain);
	get_sys_param_float(MR_GEO_GAIN_POS_X, &kpx);
	get_sys_param_float(MR_GEO_GAIN_VEL_X, &kvx);
	get_sys_param_float(MR_GEO_GAIN_POS_Y, &kpy);
	get_sys_param_float(MR_GEO_GAIN_VEL_Y, &kvy);
	get_sys_param_float(MR_GEO_GAIN_POS_Z, &kpz);
	get_sys_param_float(MR_GEO_GAIN_VEL_Z, &kvz);
	get_sys_param_float(MR_GEO_GAIN_POS_X_I, &k_tracking_i_gain[0]);
	get_sys_param_float(MR_GEO_GAIN_POS_Y_I, &k_tracking_i_gain[1]);
	get_sys_param_float(MR_GEO_GAIN_POS_Z_I, &k_tracking_i_gain[2]);
	get_sys_param_float(MR_GEO_UAV_MASS, &uav_mass);
	get_sys_param_float(MR_GEO_INERTIA_JXX, &mat_data(J)[0*3 + 0]);
	get_sys_param_float(MR_GEO_INERTIA_JYY, &mat_data(J)[1*3 + 1]);
	get_sys_param_float(MR_GEO_INERTIA_JZZ, &mat_data(J)[2*3 + 2]);
	get_sys_param_float(PWM_TO_THRUST_C1, &coeff_cmd_to_thrust[0]);
	get_sys_param_float(PWM_TO_THRUST_C2, &coeff_cmd_to_thrust[1]);
	get_sys_param_float(PWM_TO_THRUST_C3, &coeff_cmd_to_thrust[2]);
	get_sys_param_float(PWM_TO_THRUST_C4, &coeff_cmd_to_thrust[3]);
	get_sys_param_float(PWM_TO_THRUST_C5, &coeff_cmd_to_thrust[4]);
	get_sys_param_float(PWM_TO_THRUST_C6, &coeff_cmd_to_thrust[5]);
	get_sys_param_float(THRUST_TO_PWM_C1, &coeff_thrust_to_cmd[0]);
	get_sys_param_float(THRUST_TO_PWM_C2, &coeff_thrust_to_cmd[1]);
	get_sys_param_float(THRUST_TO_PWM_C3, &coeff_thrust_to_cmd[2]);
	get_sys_param_float(THRUST_TO_PWM_C4, &coeff_thrust_to_cmd[3]);
	get_sys_param_float(THRUST_TO_PWM_C5, &coeff_thrust_to_cmd[4]);
	get_sys_param_float(THRUST_TO_PWM_C6, &coeff_thrust_to_cmd[5]);
	get_sys_param_float(THRUST_MAX, &motor_thrust_max);

	set_motor_max_thrust(motor_thrust_max);
	set_motor_cmd_to_thrust_coeff(coeff_cmd_to_thrust[0], coeff_cmd_to_thrust[1], coeff_cmd_to_thrust[2],
	                              coeff_cmd_to_thrust[3], coeff_cmd_to_thrust[4], coeff_cmd_to_thrust[5]);
	set_motor_thrust_to_cmd_coeff(coeff_thrust_to_cmd[0], coeff_thrust_to_cmd[1], coeff_thrust_to_cmd[2],
	                              coeff_thrust_to_cmd[3], coeff_thrust_to_cmd[4], coeff_thrust_to_cmd[5]);
}

void estimate_uav_dynamics(float *gyro, float *moments, float *m_rot_frame)
{
	static float angular_vel_last[3] = {0.0f};
	float angular_accel[3];

	angular_accel[0] = (gyro[0] - angular_vel_last[0]) / dt;
	angular_accel[1] = (gyro[1] - angular_vel_last[1]) / dt;
	angular_accel[2] = (gyro[2] - angular_vel_last[2]) / dt;
	angular_vel_last[0] = gyro[0];
	angular_vel_last[1] = gyro[1];
	angular_vel_last[2] = gyro[2];

	lpf_first_order(angular_accel[0], &mat_data(W_dot)[0], 0.01);
	lpf_first_order(angular_accel[1], &mat_data(W_dot)[1], 0.01);
	lpf_first_order(angular_accel[2], &mat_data(W_dot)[2], 0.01);

	//J* W_dot
	MAT_MULT(&J, &W_dot, &JWdot);
	//W x JW
	MAT_MULT(&J, &W, &JW);
	cross_product_3x1(mat_data(W), mat_data(JW), mat_data(WJW));
	//M = J * W_dot + W X (J * W)
	MAT_ADD(&JWdot, &WJW, &M);

	m_rot_frame[0] = mat_data(JWdot)[0];
	m_rot_frame[1] = mat_data(JWdot)[1];
	m_rot_frame[2] = mat_data(JWdot)[2];
	moments[0] = mat_data(M)[0];
	moments[1] = mat_data(M)[1];
	moments[2] = mat_data(M)[2];
}

void reset_geometry_tracking_error_integral(void)
{
	tracking_error_integral[0] = 0.0f;
	tracking_error_integral[1] = 0.0f;
	tracking_error_integral[2] = 0.0f;
}

void geometry_manual_ctrl(euler_t *rc, float *attitude_q, float *gyro, float *output_moments,
                          bool heading_present)
{
	/* convert radio command (euler angle) to rotation matrix */
	euler_to_rotation_matrix(rc, mat_data(Rd), mat_data(Rtd));

	/* W (angular velocity) */
	mat_data(W)[0] = gyro[0];
	mat_data(W)[1] = gyro[1];
	mat_data(W)[2] = gyro[2];

	/* set Wd and Wd_dot to 0 since there is no predefined trajectory */
	mat_data(Wd)[0] = 0.0f;
	mat_data(Wd)[1] = 0.0f;
	mat_data(Wd)[2] = 0.0f;
	mat_data(Wd_dot)[0] = 0.0f;
	mat_data(Wd_dot)[1] = 0.0f;
	mat_data(Wd_dot)[2] = 0.0f;

	float _krz, _kwz; //switch between full heading control and yaw rate control

	/* switch to yaw rate control mode if no heading information provided */
	if(heading_present == false) {
		/* yaw rate control only */
		_krz = 0.0f;
		_kwz = yaw_rate_ctrl_gain;
		mat_data(Wd)[2] = rc->yaw; //set yaw rate desired value
	} else {
		_krz = krz;
		_kwz = kwz;
	}

	/* calculate attitude error eR */
	MAT_MULT(&Rtd, &R, &RtdR);
	MAT_MULT(&Rt, &Rd, &RtRd);
	MAT_SUB(&RtdR, &RtRd, &eR_mat);
	vee_map_3x3(mat_data(eR_mat), mat_data(eR));
	mat_data(eR)[0] *= 0.5f;
	mat_data(eR)[1] *= 0.5f;
	mat_data(eR)[2] *= 0.5f;

	/* calculate attitude rate error eW */
	//MAT_MULT(&Rt, &Rd, &RtRd); //the term is duplicated
	MAT_MULT(&RtRd, &Wd, &RtRdWd);
	MAT_SUB(&W, &RtRdWd, &eW);

	/* calculate the inertia feedfoward term */
	//W x JW
	MAT_MULT(&J, &W, &JW);
	cross_product_3x1(mat_data(W), mat_data(JW), mat_data(WJW));
	mat_data(inertia_effect)[0] = mat_data(WJW)[0];
	mat_data(inertia_effect)[1] = mat_data(WJW)[1];
	mat_data(inertia_effect)[2] = mat_data(WJW)[2];

#if 0   /* inertia feedfoward term for motion planning (trajectory is known) */
	/* calculate inertia effect (trajectory is defined, Wd and Wd_dot are not zero) */
	//W * R^T * Rd * Wd
	hat_map_3x3(mat_data(W), mat_data(W_hat));
	MAT_MULT(&W_hat, &Rt, &WRt);
	MAT_MULT(&WRt, &Rd, &WRtRd);
	MAT_MULT(&WRtRd, &Wd, &WRtRdWd);
	//R^T * Rd * Wd_dot
	//MAT_MULT(&Rt, &Rd, &RtRd); //the term is duplicated
	MAT_MULT(&RtRd, &Wd_dot, &RtRdWddot);
	//(W * R^T * Rd * Wd) - (R^T * Rd * Wd_dot)
	MAT_SUB(&WRtRdWd, &RtRdWddot, &WRtRdWd_RtRdWddot);
	//J*[(W * R^T * Rd * Wd) - (R^T * Rd * Wd_dot)]
	MAT_MULT(&J, &WRtRdWd_RtRdWddot, &J_WRtRdWd_RtRdWddot);
	//inertia effect = (W x JW) - J*[(W * R^T * Rd * Wd) - (R^T * Rd * Wd_dot)]
	MAT_SUB(&WJW, &J_WRtRdWd_RtRdWddot, &inertia_effect);

#endif

	/* control input M1, M2, M3 */
	output_moments[0] = -krx*mat_data(eR)[0] -kwx*mat_data(eW)[0] + mat_data(inertia_effect)[0];
	output_moments[1] = -kry*mat_data(eR)[1] -kwy*mat_data(eW)[1] + mat_data(inertia_effect)[1];
	output_moments[2] = -_krz*mat_data(eR)[2] -_kwz*mat_data(eW)[2] + mat_data(inertia_effect)[2];
}

void geometry_tracking_ctrl(euler_t *rc, float *attitude_q, float *gyro,
                            float *pos_des_enu, float *vel_des_enu, float *accel_ff_enu,
                            float *curr_pos_ned, float *curr_vel_ned, float *output_moments,
                            float *output_force, bool manual_flight)
{
	/* ex = x - xd */
	float pos_des_ned[3];
	assign_vector_3x1_enu_to_ned(pos_des_ned, pos_des_enu);
	pos_error[0] = curr_pos_ned[0] - pos_des_ned[0];
	pos_error[1] = curr_pos_ned[1] - pos_des_ned[1];
	pos_error[2] = curr_pos_ned[2] - pos_des_ned[2];

	/* ev = v - vd */
	float vel_des_ned[3];
	assign_vector_3x1_enu_to_ned(vel_des_ned, vel_des_enu);
	vel_error[0] = curr_vel_ned[0] - vel_des_ned[0];
	vel_error[1] = curr_vel_ned[1] - vel_des_ned[1];
	vel_error[2] = curr_vel_ned[2] - vel_des_ned[2];

	float force_ff_ned[3];
	float accel_ff_ned[3];
	assign_vector_3x1_enu_to_ned(accel_ff_ned, accel_ff_enu);
	force_ff_ned[0] = uav_mass * accel_ff_ned[0];
	force_ff_ned[1] = uav_mass * accel_ff_ned[1];
	force_ff_ned[2] = uav_mass * accel_ff_ned[2];

	tracking_error_integral[0] += k_tracking_i_gain[0] * (pos_error[0]) * dt;
	tracking_error_integral[1] += k_tracking_i_gain[1] * (pos_error[1]) * dt;
	tracking_error_integral[2] += k_tracking_i_gain[2] * (pos_error[2]) * dt;

	bound_float(&tracking_error_integral[0], 150, -150);
	bound_float(&tracking_error_integral[1], 150, -150);
	bound_float(&tracking_error_integral[2], 50, -50);

	mat_data(kxex_kvev_mge3_mxd_dot_dot)[0] = -kpx*pos_error[0] - kvx*vel_error[0] +
	                force_ff_ned[0] - tracking_error_integral[0];
	mat_data(kxex_kvev_mge3_mxd_dot_dot)[1] = -kpy*pos_error[1] - kvy*vel_error[1] +
	                force_ff_ned[1] - tracking_error_integral[1];
	mat_data(kxex_kvev_mge3_mxd_dot_dot)[2] = -kpz*pos_error[2] - kvz*vel_error[2] +
	                force_ff_ned[2] - tracking_error_integral[2] -
	                uav_mass * 9.81;

	/* calculate the denominator of b3d */
	float b3d_denominator; //caution: this term should not be 0
	norm_3x1(mat_data(kxex_kvev_mge3_mxd_dot_dot), &b3d_denominator);
	b3d_denominator = -1.0f / b3d_denominator;

	if(manual_flight == true) {
		/* enable altitude control only, control roll and pitch manually */
		//convert radio command (euler angle) to rotation matrix
		euler_to_rotation_matrix(rc, mat_data(Rd), mat_data(Rtd));
	} else {
		/* enable tracking control for x and y axis */
		//b1d
		mat_data(b1d)[0] = arm_cos_f32(rc->yaw);
		mat_data(b1d)[1] = arm_sin_f32(rc->yaw);
		mat_data(b1d)[2] = 0.0f;
		//b3d = -kxex_kvev_mge3_mxd_dot_dot / ||kxex_kvev_mge3_mxd_dot_dot||
		mat_data(b3d)[0] = mat_data(kxex_kvev_mge3_mxd_dot_dot)[0] * b3d_denominator;
		mat_data(b3d)[1] = mat_data(kxex_kvev_mge3_mxd_dot_dot)[1] * b3d_denominator;
		mat_data(b3d)[2] = mat_data(kxex_kvev_mge3_mxd_dot_dot)[2] * b3d_denominator;
		//b2d = b3d X b1d / ||b3d X b1d||
		cross_product_3x1(mat_data(b3d), mat_data(b1d), mat_data(b2d));
		normalize_3x1(mat_data(b2d));
		/* proj[b1d] = b2d X b3d */
		cross_product_3x1(mat_data(b2d), mat_data(b3d), mat_data(b1d));

		//Rd = [b1d; b3d X b1d; b3d]
		mat_data(Rd)[0*3 + 0] = mat_data(b1d)[0];
		mat_data(Rd)[1*3 + 0] = mat_data(b1d)[1];
		mat_data(Rd)[2*3 + 0] = mat_data(b1d)[2];
		mat_data(Rd)[0*3 + 1] = mat_data(b2d)[0];
		mat_data(Rd)[1*3 + 1] = mat_data(b2d)[1];
		mat_data(Rd)[2*3 + 1] = mat_data(b2d)[2];
		mat_data(Rd)[0*3 + 2] = mat_data(b3d)[0];
		mat_data(Rd)[1*3 + 2] = mat_data(b3d)[1];
		mat_data(Rd)[2*3 + 2] = mat_data(b3d)[2];

		//transpose(Rd)
		mat_data(Rtd)[0*3 + 0] = mat_data(Rd)[0*3 + 0];
		mat_data(Rtd)[1*3 + 0] = mat_data(Rd)[0*3 + 1];
		mat_data(Rtd)[2*3 + 0] = mat_data(Rd)[0*3 + 2];
		mat_data(Rtd)[0*3 + 1] = mat_data(Rd)[1*3 + 0];
		mat_data(Rtd)[1*3 + 1] = mat_data(Rd)[1*3 + 1];
		mat_data(Rtd)[2*3 + 1] = mat_data(Rd)[1*3 + 2];
		mat_data(Rtd)[0*3 + 2] = mat_data(Rd)[2*3 + 0];
		mat_data(Rtd)[1*3 + 2] = mat_data(Rd)[2*3 + 1];
		mat_data(Rtd)[2*3 + 2] = mat_data(Rd)[2*3 + 2];
	}

	/* R * e3 */
	mat_data(Re3)[0] = mat_data(R)[0*3 + 2];
	mat_data(Re3)[1] = mat_data(R)[1*3 + 2];
	mat_data(Re3)[2] = mat_data(R)[2*3 + 2];
	/* f = -(-kx * ex - kv * ev - mge3 + m * x_d_dot_dot) . (R * e3) */
	float neg_kxex_kvev_mge3_mxd_dot_dot[3];
	neg_kxex_kvev_mge3_mxd_dot_dot[0] = -mat_data(kxex_kvev_mge3_mxd_dot_dot)[0];
	neg_kxex_kvev_mge3_mxd_dot_dot[1] = -mat_data(kxex_kvev_mge3_mxd_dot_dot)[1];
	neg_kxex_kvev_mge3_mxd_dot_dot[2] = -mat_data(kxex_kvev_mge3_mxd_dot_dot)[2];
	arm_dot_prod_f32(neg_kxex_kvev_mge3_mxd_dot_dot, mat_data(Re3), 3, output_force);

	/* W (angular velocity) */
	mat_data(W)[0] = gyro[0];
	mat_data(W)[1] = gyro[1];
	mat_data(W)[2] = gyro[2];

	/* set Wd and Wd_dot to 0 since there is no predefined trajectory */
	mat_data(Wd)[0] = 0.0f;
	mat_data(Wd)[1] = 0.0f;
	mat_data(Wd)[2] = 0.0f;
	mat_data(Wd_dot)[0] = 0.0f;
	mat_data(Wd_dot)[1] = 0.0f;
	mat_data(Wd_dot)[2] = 0.0f;

	/* calculate attitude error eR */
	MAT_MULT(&Rtd, &R, &RtdR);
	MAT_MULT(&Rt, &Rd, &RtRd);
	MAT_SUB(&RtdR, &RtRd, &eR_mat);
	vee_map_3x3(mat_data(eR_mat), mat_data(eR));
	mat_data(eR)[0] *= 0.5f;
	mat_data(eR)[1] *= 0.5f;
	mat_data(eR)[2] *= 0.5f;

	/* calculate attitude rate error eW */
	//MAT_MULT(&Rt, &Rd, &RtRd); //the term is duplicated
	MAT_MULT(&RtRd, &Wd, &RtRdWd);
	MAT_SUB(&W, &RtRdWd, &eW);

	/* calculate the inertia feedfoward term */
	//W x JW
	MAT_MULT(&J, &W, &JW);
	cross_product_3x1(mat_data(W), mat_data(JW), mat_data(WJW));
	mat_data(inertia_effect)[0] = mat_data(WJW)[0];
	mat_data(inertia_effect)[1] = mat_data(WJW)[1];
	mat_data(inertia_effect)[2] = mat_data(WJW)[2];

#if (SELECT_CONTROLLER_ESTIMATOR == CONTROLLER_ESTIMATION_USE_ORIGIN)
#pragma "ORIGIN"
	/* control input M1, M2, M3 */
	output_moments[0] = -krx*mat_data(eR)[0] -kwx*mat_data(eW)[0] + mat_data(inertia_effect)[0];
	output_moments[1] = -kry*mat_data(eR)[1] -kwy*mat_data(eW)[1] + mat_data(inertia_effect)[1];
	output_moments[2] = -krz*mat_data(eR)[2] -kwz*mat_data(eW)[2] + mat_data(inertia_effect)[2];
#elif (SELECT_CONTROLLER_ESTIMATOR == CONTROLLER_ESTIMATION_USE_NONE)
	
#pragma "NONE"
	/*adaptive*/	
	//Y1 of Y
	// [0	-f	]
	// [f	0	]
	// [0	0	]
	mat_data(Y)[0*8 + 0] = 0; 
	mat_data(Y)[1*8 + 0] = *output_force; 
	mat_data(Y)[2*8 + 0] = 0;
	mat_data(Y)[0*8 + 1] = -*output_force; 
	mat_data(Y)[1*8 + 1] = 0; 
	mat_data(Y)[2*8 + 1] = 0;
	//Y2 of Y
	//wb: omega_bar
	//w: angular velocity
	//[-wb1		,-w1*w2	,w1*w2	,-wb1-w0*w2	,-wb2+w0*w1	,-w2^2+w1^2	]
	//[w0*w2	,-wb1	,-w0*w2	,-wb0+w1*w2	,-w0^2+w2^2	,-wb2-w0*w1	]
	//[-w0*w1	,w0*w1	,--wb2	,-w1^2+w0^2	,-wb0-w1*w2	,-wb1+w0*w2	]
	mat_data(Y)[0*8 + 2] = 0; 
	mat_data(Y)[1*8 + 2] = mat_data(W)[0]*mat_data(W)[2];
	mat_data(Y)[2*8 + 2] = -mat_data(W)[0]*mat_data(W)[1];
	mat_data(Y)[0*8 + 3] = -mat_data(W)[1]*mat_data(W)[2];
	mat_data(Y)[1*8 + 3] = 0;
	mat_data(Y)[2*8 + 3] = mat_data(W)[0]*mat_data(W)[1];
	mat_data(Y)[0*8 + 4] = mat_data(W)[1]*mat_data(W)[2];
	mat_data(Y)[1*8 + 4] = -mat_data(W)[0]*mat_data(W)[2];
	mat_data(Y)[2*8 + 4] = 0;
	mat_data(Y)[0*8 + 5] = -mat_data(W)[0]*mat_data(W)[2];
	mat_data(Y)[1*8 + 5] = mat_data(W)[1]*mat_data(W)[2];
	mat_data(Y)[2*8 + 5] = -mat_data(W)[1]*mat_data(W)[1] +mat_data(W)[0]*mat_data(W)[0];
	mat_data(Y)[0*8 + 6] = mat_data(W)[0]*mat_data(W)[1];
	mat_data(Y)[1*8 + 6] = -mat_data(W)[0]*mat_data(W)[0] +mat_data(W)[2]*mat_data(W)[2];
	mat_data(Y)[2*8 + 6] = -mat_data(W)[1]*mat_data(W)[2];
	mat_data(Y)[0*8 + 7] = -mat_data(W)[2]*mat_data(W)[2] +mat_data(W)[1]*mat_data(W)[1];
	mat_data(Y)[1*8 + 7] = -mat_data(W)[0]*mat_data(W)[1];
	mat_data(Y)[2*8 + 7] = mat_data(W)[0]*mat_data(W)[2];

	//Y_theta
	MAT_MULT(&Y, &theta, &Ytheta);
	
	/* control input M1, M2, M3 */
	output_moments[0] = -krx*mat_data(eR)[0] -kwx*mat_data(eW)[0] +mat_data(Ytheta)[0]; 
	output_moments[1] = -kry*mat_data(eR)[1] -kwy*mat_data(eW)[1] +mat_data(Ytheta)[1];
	output_moments[2] = -krz*mat_data(eR)[2] -kwz*mat_data(eW)[2] +mat_data(Ytheta)[2];
#elif (SELECT_CONTROLLER_ESTIMATOR == CONTROLLER_ESTIMATION_USE_ADAPTIVE)
	
	/*adaptive*/	
	//Y1 of Y
	// [0	-f	]
	// [f	0	]
	// [0	0	]
	mat_data(Y)[0*8 + 0] = 0; 
	mat_data(Y)[1*8 + 0] = *output_force; 
	mat_data(Y)[2*8 + 0] = 0;
	mat_data(Y)[0*8 + 1] = -*output_force; 
	mat_data(Y)[1*8 + 1] = 0; 
	mat_data(Y)[2*8 + 1] = 0;
	//Y2 of Y
	//wb: omega_bar
	//w: angular velocity
	//[-wb1		,-w1*w2	,w1*w2	,-wb1-w0*w2	,-wb2+w0*w1	,-w2^2+w1^2	]
	//[w0*w2	,-wb1	,-w0*w2	,-wb0+w1*w2	,-w0^2+w2^2	,-wb2-w0*w1	]
	//[-w0*w1	,w0*w1	,--wb2	,-w1^2+w0^2	,-wb0-w1*w2	,-wb1+w0*w2	]
	mat_data(Y)[0*8 + 2] = 0; 
	mat_data(Y)[1*8 + 2] = mat_data(W)[0]*mat_data(W)[2];
	mat_data(Y)[2*8 + 2] = -mat_data(W)[0]*mat_data(W)[1];
	mat_data(Y)[0*8 + 3] = -mat_data(W)[1]*mat_data(W)[2];
	mat_data(Y)[1*8 + 3] = 0;
	mat_data(Y)[2*8 + 3] = mat_data(W)[0]*mat_data(W)[1];
	mat_data(Y)[0*8 + 4] = mat_data(W)[1]*mat_data(W)[2];
	mat_data(Y)[1*8 + 4] = -mat_data(W)[0]*mat_data(W)[2];
	mat_data(Y)[2*8 + 4] = 0;
	mat_data(Y)[0*8 + 5] = -mat_data(W)[0]*mat_data(W)[2];
	mat_data(Y)[1*8 + 5] = mat_data(W)[1]*mat_data(W)[2];
	mat_data(Y)[2*8 + 5] = -mat_data(W)[1]*mat_data(W)[1] +mat_data(W)[0]*mat_data(W)[0];
	mat_data(Y)[0*8 + 6] = mat_data(W)[0]*mat_data(W)[1];
	mat_data(Y)[1*8 + 6] = -mat_data(W)[0]*mat_data(W)[0] +mat_data(W)[2]*mat_data(W)[2];
	mat_data(Y)[2*8 + 6] = -mat_data(W)[1]*mat_data(W)[2];
	mat_data(Y)[0*8 + 7] = -mat_data(W)[2]*mat_data(W)[2] +mat_data(W)[1]*mat_data(W)[1];
	mat_data(Y)[1*8 + 7] = -mat_data(W)[0]*mat_data(W)[1];
	mat_data(Y)[2*8 + 7] = mat_data(W)[0]*mat_data(W)[2];
	//Yt*(eW_c2eR)	
	MAT_TRANS(&Y, &Yt);
	mat_data(eW_c2eR)[0] = c2*mat_data(eR)[0] + mat_data(eW)[0];
	mat_data(eW_c2eR)[1] = c2*mat_data(eR)[1] + mat_data(eW)[1];
	mat_data(eW_c2eR)[2] = c2*mat_data(eR)[2] + mat_data(eW)[2];	
	MAT_MULT(&Yt,&eW_c2eR,&Yt_eW_c2eR);
	
	//adaptivetheta_hat_dot 
	for(int i = 0 ; i<8 ; i++){
		mat_data(theta_hat_dot)[i] = -adaptive_gamma[i]*mat_data(Yt_eW_c2eR)[i];
	}
	//theta update
	mat_data(theta)[0] = mat_data(theta)[0] + mat_data(theta_hat_dot)[0];
	mat_data(theta)[1] = mat_data(theta)[1] + mat_data(theta_hat_dot)[1];
	mat_data(theta)[2] = mat_data(theta)[2] + mat_data(theta_hat_dot)[2];
	mat_data(theta)[3] = mat_data(theta)[3] + mat_data(theta_hat_dot)[3];
	mat_data(theta)[4] = mat_data(theta)[4] + mat_data(theta_hat_dot)[4];
	mat_data(theta)[5] = mat_data(theta)[5] + mat_data(theta_hat_dot)[5];
	mat_data(theta)[6] = mat_data(theta)[6] + mat_data(theta_hat_dot)[6];
	mat_data(theta)[7] = mat_data(theta)[7] + mat_data(theta_hat_dot)[7];

	//Y_theta
	MAT_MULT(&Y, &theta, &Ytheta);
	
	/* control input M1, M2, M3 */
	output_moments[0] = -krx*mat_data(eR)[0] -kwx*mat_data(eW)[0] +mat_data(Ytheta)[0]; 
	output_moments[1] = -kry*mat_data(eR)[1] -kwy*mat_data(eW)[1] +mat_data(Ytheta)[1];
	output_moments[2] = -krz*mat_data(eR)[2] -kwz*mat_data(eW)[2] +mat_data(Ytheta)[2];
#elif (SELECT_CONTROLLER_ESTIMATOR == CONTROLLER_ESTIMATION_USE_ICL)
#pragma "ICL"	
	/*adaptive*/	
	//Y1 of Y
	// [0	-f	]
	// [f	0	]
	// [0	0	]
	//
	mat_data(Y)[0*8 + 0] = 0; 
	mat_data(Y)[1*8 + 0] = *output_force; 
	mat_data(Y)[2*8 + 0] = 0;
	mat_data(Y)[0*8 + 1] = -*output_force; 
	mat_data(Y)[1*8 + 1] = 0; 
	mat_data(Y)[2*8 + 1] = 0;
	//Y2 of Y
	//wb: omega_bar
	//w: angular velocity
	//[-wb1		,-w1*w2	,w1*w2	,-wb1-w0*w2	,-wb2+w0*w1	,-w2^2+w1^2	]
	//[w0*w2	,-wb1	,-w0*w2	,-wb0+w1*w2	,-w0^2+w2^2	,-wb2-w0*w1	]
	//[-w0*w1	,w0*w1	,--wb2	,-w1^2+w0^2	,-wb0-w1*w2	,-wb1+w0*w2	]
	mat_data(Y)[0*8 + 2] = 0; 
	mat_data(Y)[1*8 + 2] = mat_data(W)[0]*mat_data(W)[2];
	mat_data(Y)[2*8 + 2] = -mat_data(W)[0]*mat_data(W)[1];
	mat_data(Y)[0*8 + 3] = -mat_data(W)[1]*mat_data(W)[2];
	mat_data(Y)[1*8 + 3] = 0;
	mat_data(Y)[2*8 + 3] = mat_data(W)[0]*mat_data(W)[1];
	mat_data(Y)[0*8 + 4] = mat_data(W)[1]*mat_data(W)[2];
	mat_data(Y)[1*8 + 4] = -mat_data(W)[0]*mat_data(W)[2];
	mat_data(Y)[2*8 + 4] = 0;
	mat_data(Y)[0*8 + 5] = -mat_data(W)[0]*mat_data(W)[2];
	mat_data(Y)[1*8 + 5] = mat_data(W)[1]*mat_data(W)[2];
	mat_data(Y)[2*8 + 5] = -mat_data(W)[1]*mat_data(W)[1] +mat_data(W)[0]*mat_data(W)[0];
	mat_data(Y)[0*8 + 6] = mat_data(W)[0]*mat_data(W)[1];
	mat_data(Y)[1*8 + 6] = -mat_data(W)[0]*mat_data(W)[0] +mat_data(W)[2]*mat_data(W)[2];
	mat_data(Y)[2*8 + 6] = -mat_data(W)[1]*mat_data(W)[2];
	mat_data(Y)[0*8 + 7] = -mat_data(W)[2]*mat_data(W)[2] +mat_data(W)[1]*mat_data(W)[1];
	mat_data(Y)[1*8 + 7] = -mat_data(W)[0]*mat_data(W)[1];
	mat_data(Y)[2*8 + 7] = mat_data(W)[0]*mat_data(W)[2];
	// update adaptive theta
	//Yt*(eW_c2eR)	
	MAT_TRANS(&Y, &Yt);
	mat_data(eW_c2eR)[0] = c2*mat_data(eR)[0] + mat_data(eW)[0];
	mat_data(eW_c2eR)[1] = c2*mat_data(eR)[1] + mat_data(eW)[1];
	mat_data(eW_c2eR)[2] = c2*mat_data(eR)[2] + mat_data(eW)[2];	
	MAT_MULT(&Yt,&eW_c2eR,&Yt_eW_c2eR);
	
	//theta_hat_dot 
	for(int i = 0 ; i<8 ; i++){
		mat_data(adaptive_theta_hat_dot)[i] = -adaptive_gamma[i]*mat_data(Yt_eW_c2eR)[i];
	}

	/*ICL*/
	//Y1 of Y in icl	
	// [0	-f	]
	// [f	0	]
	// [0	0	]
	//float ave_force = (output_force_last + output_force_last2 + output_force_last3 + output_force_last4 + output_force_last5 )/5;
	mat_data(Y1_Y_W)[0*8 + 0] = 0; 
	mat_data(Y1_Y_W)[1*8 + 0] = output_force_last5; 
	mat_data(Y1_Y_W)[2*8 + 0] = 0;
	mat_data(Y1_Y_W)[0*8 + 1] = -output_force_last5; 
	mat_data(Y1_Y_W)[1*8 + 1] = 0; 
	mat_data(Y1_Y_W)[2*8 + 1] = 0;
	//Y W of Y in icl	
	//w: angular velocity
	//[0		,-w1*w2	,w1*w2	,-w0*w2		,w0*w1		,-w2^2+w1^2	]
	//[w0*w2	,0		,-w0*w2	,w1*w2		,-w0^2+w2^2	,-w0*w1		]
	//[-w0*w1	,w0*w1	,0		,-w1^2+w0^2	,-w1*w2		,w0*w2		]
	/*
	float W0 = (mat_data(W_last)[0] + mat_data(W_last2)[0] + mat_data(W_last3)[0] + mat_data(W_last4)[0] + mat_data(W_last5)[0])/5;
	float W1 = (mat_data(W_last)[1] + mat_data(W_last2)[1] + mat_data(W_last3)[1] + mat_data(W_last4)[1] + mat_data(W_last5)[1])/5;
	float W2 = (mat_data(W_last)[2] + mat_data(W_last2)[2] + mat_data(W_last3)[2] + mat_data(W_last4)[2] + mat_data(W_last5)[2])/5;
	
	mat_data(Y1_Y_W)[0*8 + 2] = 0; 
	mat_data(Y1_Y_W)[1*8 + 2] = W0*W2;
	mat_data(Y1_Y_W)[2*8 + 2] = -W0*W1;
	mat_data(Y1_Y_W)[0*8 + 3] = -W1*W2;
	mat_data(Y1_Y_W)[1*8 + 3] = 0;
	mat_data(Y1_Y_W)[2*8 + 3] = W0*W1;
	mat_data(Y1_Y_W)[0*8 + 4] = W1*W2;
	mat_data(Y1_Y_W)[1*8 + 4] = -W0*W2;
	mat_data(Y1_Y_W)[2*8 + 4] = 0;
	mat_data(Y1_Y_W)[0*8 + 5] = -W0*W2;
	mat_data(Y1_Y_W)[1*8 + 5] = W1*W2;
	mat_data(Y1_Y_W)[2*8 + 5] = -W1*W1 + W0*W0;
	mat_data(Y1_Y_W)[0*8 + 6] = W0*W1;
	mat_data(Y1_Y_W)[1*8 + 6] = -W0*W0 + W2*W2;
	mat_data(Y1_Y_W)[2*8 + 6] = -W1*W2;
	mat_data(Y1_Y_W)[0*8 + 7] = -W2*W2 + W1*W1;
	mat_data(Y1_Y_W)[1*8 + 7] = -W0*W1;
	mat_data(Y1_Y_W)[2*8 + 7] = W0*W2;
	*/
	mat_data(Y1_Y_W)[0*8 + 2] = 0; 
	mat_data(Y1_Y_W)[1*8 + 2] = mat_data(W)[0]*mat_data(W)[2];
	mat_data(Y1_Y_W)[2*8 + 2] = -mat_data(W)[0]*mat_data(W)[1];
	mat_data(Y1_Y_W)[0*8 + 3] = -mat_data(W)[1]*mat_data(W)[2];
	mat_data(Y1_Y_W)[1*8 + 3] = 0;
	mat_data(Y1_Y_W)[2*8 + 3] = mat_data(W)[0]*mat_data(W)[1];
	mat_data(Y1_Y_W)[0*8 + 4] = mat_data(W)[1]*mat_data(W)[2];
	mat_data(Y1_Y_W)[1*8 + 4] = -mat_data(W)[0]*mat_data(W)[2];
	mat_data(Y1_Y_W)[2*8 + 4] = 0;
	mat_data(Y1_Y_W)[0*8 + 5] = -mat_data(W)[0]*mat_data(W)[2];
	mat_data(Y1_Y_W)[1*8 + 5] = mat_data(W)[1]*mat_data(W)[2];
	mat_data(Y1_Y_W)[2*8 + 5] = -mat_data(W)[1]*mat_data(W)[1] +mat_data(W)[0]*mat_data(W)[0];
	mat_data(Y1_Y_W)[0*8 + 6] = mat_data(W)[0]*mat_data(W)[1];
	mat_data(Y1_Y_W)[1*8 + 6] = -mat_data(W)[0]*mat_data(W)[0] +mat_data(W)[2]*mat_data(W)[2];
	mat_data(Y1_Y_W)[2*8 + 6] = -mat_data(W)[1]*mat_data(W)[2];
	mat_data(Y1_Y_W)[0*8 + 7] = -mat_data(W)[2]*mat_data(W)[2] +mat_data(W)[1]*mat_data(W)[1];
	mat_data(Y1_Y_W)[1*8 + 7] = -mat_data(W)[0]*mat_data(W)[1];
	mat_data(Y1_Y_W)[2*8 + 7] = mat_data(W)[0]*mat_data(W)[2];


#if 1
	/* generate Inertia filter*/
	MAT_ALLOC(J_gt, 3, 3);
	MAT_INIT(J_gt, 3, 3);
	MAT_ALLOC(J_gt_inv, 3, 3);
	MAT_INIT(J_gt_inv, 3, 3);
	MAT_ALLOC(f_cog, 3, 1);
	MAT_INIT(f_cog, 3, 1);
	/* define J_gt */
#if 1
	mat_data(J_gt)[0*3 + 0] = mat_data(theta)[2]; 
	mat_data(J_gt)[1*3 + 0] = mat_data(theta)[5]; 
	mat_data(J_gt)[2*3 + 0] = mat_data(theta)[6];
	mat_data(J_gt)[0*3 + 1] = mat_data(theta)[5]; 
	mat_data(J_gt)[1*3 + 1] = mat_data(theta)[3]; 
	mat_data(J_gt)[2*3 + 1] = mat_data(theta)[7];
	mat_data(J_gt)[0*3 + 2] = mat_data(theta)[6]; 
	mat_data(J_gt)[1*3 + 2] = mat_data(theta)[7]; 
	mat_data(J_gt)[2*3 + 2] = mat_data(theta)[4];
#else
	mat_data(J_gt)[0*3 + 0] = 0.0133; 
	mat_data(J_gt)[1*3 + 0] = 0.0021; 
	mat_data(J_gt)[2*3 + 0] = 0;
	mat_data(J_gt)[0*3 + 1] = 0.0021; 
	mat_data(J_gt)[1*3 + 1] = 0.0134; 
	mat_data(J_gt)[2*3 + 1] = 0;
	mat_data(J_gt)[0*3 + 2] = 0; 
	mat_data(J_gt)[1*3 + 2] = 0; 
	mat_data(J_gt)[2*3 + 2] = 0.0228;
#endif
	/*define J_gt_inv */
	MAT_INV(&J_gt, &J_gt_inv);
	/*define r_cog */
	mat_data(f_cog)[0] = -mat_data(theta)[1]*output_force_last5; 
	mat_data(f_cog)[1] = mat_data(theta)[0]*output_force_last5; 
	mat_data(f_cog)[2] = 0;
	/*genetate WJW */	
	MAT_MULT(&J_gt, &W, &JW);
	cross_product_3x1(mat_data(W), mat_data(JW), mat_data(WJW));
	/* generate W_dot_hat */
	MAT_ALLOC(W_dot_hat, 3, 1);
	MAT_INIT(W_dot_hat, 3, 1);
	MAT_ALLOC(m_f_cog, 3, 1);
	MAT_INIT(m_f_cog, 3, 1);
	MAT_ALLOC(m_f_cog_WJW, 3, 1);
	MAT_INIT(m_f_cog_WJW, 3, 1);

	MAT_SUB(&M_last5, &f_cog, &m_f_cog);
	MAT_SUB(&m_f_cog, &WJW, &m_f_cog_WJW);
	MAT_MULT(&J_gt_inv, &m_f_cog_WJW, &W_dot_hat);
	/* generate W_filt */
	MAT_ALLOC(W_filt, 3, 1);
	MAT_INIT(W_filt, 3, 1);
	float k_omega = 0.8;
	mat_data(W_filt)[0] = (1.0f - k_omega)*mat_data(W)[0] + k_omega*(mat_data(W_last5)[0] + mat_data(W_dot_hat)[0]*5*dt);
	mat_data(W_filt)[1] = (1.0f - k_omega)*mat_data(W)[1] + k_omega*(mat_data(W_last5)[1] + mat_data(W_dot_hat)[1]*5*dt);
	mat_data(W_filt)[2] = (1.0f - k_omega)*mat_data(W)[2] + k_omega*(mat_data(W_last5)[2] + mat_data(W_dot_hat)[2]*5*dt);
	/*common ICL regression*/
	MAT_SUB(&W_filt, &W_last5, &W_dot);
	mat_data(W_dot_Matrix)[0*8 + 0] = 0; 
	mat_data(W_dot_Matrix)[1*8 + 0] = 0;
	mat_data(W_dot_Matrix)[2*8 + 0] = 0;
	mat_data(W_dot_Matrix)[0*8 + 1] = 0; 
	mat_data(W_dot_Matrix)[1*8 + 1] = 0;
	mat_data(W_dot_Matrix)[2*8 + 1] = 0;
	mat_data(W_dot_Matrix)[0*8 + 2] = mat_data(W_dot)[0]; 
	mat_data(W_dot_Matrix)[1*8 + 2] = 0;
	mat_data(W_dot_Matrix)[2*8 + 2] = 0;
	mat_data(W_dot_Matrix)[0*8 + 3] = 0;
	mat_data(W_dot_Matrix)[1*8 + 3] = mat_data(W_dot)[1];
	mat_data(W_dot_Matrix)[2*8 + 3] = 0;
	mat_data(W_dot_Matrix)[0*8 + 4] = 0;
	mat_data(W_dot_Matrix)[1*8 + 4] = 0;
	mat_data(W_dot_Matrix)[2*8 + 4] = mat_data(W_dot)[2];
	mat_data(W_dot_Matrix)[0*8 + 5] = mat_data(W_dot)[1];
	mat_data(W_dot_Matrix)[1*8 + 5] = mat_data(W_dot)[0];
	mat_data(W_dot_Matrix)[2*8 + 5] = 0;
	mat_data(W_dot_Matrix)[0*8 + 6] = mat_data(W_dot)[2];
	mat_data(W_dot_Matrix)[1*8 + 6] = 0;
	mat_data(W_dot_Matrix)[2*8 + 6] = mat_data(W_dot)[0];
	mat_data(W_dot_Matrix)[0*8 + 7] = 0;
	mat_data(W_dot_Matrix)[1*8 + 7] = mat_data(W_dot)[2];
	mat_data(W_dot_Matrix)[2*8 + 7] = mat_data(W_dot)[1];
#else
	MAT_SUB(&W, &W_last5, &W_dot);
	mat_data(W_dot_Matrix)[0*8 + 0] = 0; 
	mat_data(W_dot_Matrix)[1*8 + 0] = 0;
	mat_data(W_dot_Matrix)[2*8 + 0] = 0;
	mat_data(W_dot_Matrix)[0*8 + 1] = 0; 
	mat_data(W_dot_Matrix)[1*8 + 1] = 0;
	mat_data(W_dot_Matrix)[2*8 + 1] = 0;
	mat_data(W_dot_Matrix)[0*8 + 2] = mat_data(W_dot)[0]; 
	mat_data(W_dot_Matrix)[1*8 + 2] = 0;
	mat_data(W_dot_Matrix)[2*8 + 2] = 0;
	mat_data(W_dot_Matrix)[0*8 + 3] = 0;
	mat_data(W_dot_Matrix)[1*8 + 3] = mat_data(W_dot)[1];
	mat_data(W_dot_Matrix)[2*8 + 3] = 0;
	mat_data(W_dot_Matrix)[0*8 + 4] = 0;
	mat_data(W_dot_Matrix)[1*8 + 4] = 0;
	mat_data(W_dot_Matrix)[2*8 + 4] = mat_data(W_dot)[2];
	mat_data(W_dot_Matrix)[0*8 + 5] = mat_data(W_dot)[1];
	mat_data(W_dot_Matrix)[1*8 + 5] = mat_data(W_dot)[0];
	mat_data(W_dot_Matrix)[2*8 + 5] = 0;
	mat_data(W_dot_Matrix)[0*8 + 6] = mat_data(W_dot)[2];
	mat_data(W_dot_Matrix)[1*8 + 6] = 0;
	mat_data(W_dot_Matrix)[2*8 + 6] = mat_data(W_dot)[0];
	mat_data(W_dot_Matrix)[0*8 + 7] = 0;
	mat_data(W_dot_Matrix)[1*8 + 7] = mat_data(W_dot)[2];
	mat_data(W_dot_Matrix)[2*8 + 7] = mat_data(W_dot)[1];
#endif
	//Y1_Y_Wdt = Y1_Y
	MAT_SCALE(&Y1_Y_W, dt*5, &Y1_Y_Wdt);
	
	//update sigma array
	//y_cl & Mdt
	/*
	MAT_ALLOC(M_ave,3,1);	
	MAT_INIT(M_ave,3,1);	
	mat_data(M_ave)[0] = (mat_data(M_last)[0] + mat_data(M_last2)[0] + mat_data(M_last3)[0] + mat_data(M_last4)[0] + mat_data(M_last5)[0])/5;
	mat_data(M_ave)[1] = (mat_data(M_last)[1] + mat_data(M_last2)[1] + mat_data(M_last3)[1] + mat_data(M_last4)[1] + mat_data(M_last5)[1])/5;
	mat_data(M_ave)[2] = (mat_data(M_last)[2] + mat_data(M_last2)[2] + mat_data(M_last3)[2] + mat_data(M_last4)[2] + mat_data(M_last5)[2])/5;
	*/
	MAT_ADD(&Y1_Y_Wdt ,&W_dot_Matrix ,&sigma_array[ICL_sigma_index].y_cl );
	MAT_SCALE(&M_last5, dt*5, &sigma_array[ICL_sigma_index].Mb);
	ICL_sigma_index ++;
	ICL_sigma_index %= ICL_N;

//============= calculate ICL sigma control term ===============
	mat_data(ICL_control_term)[0] = 0;
	mat_data(ICL_control_term)[1] = 0;
	mat_data(ICL_control_term)[2] = 0;
	mat_data(ICL_control_term)[3] = 0;
	mat_data(ICL_control_term)[4] = 0;
	mat_data(ICL_control_term)[5] = 0;
	mat_data(ICL_control_term)[6] = 0;
	mat_data(ICL_control_term)[7] = 0;

	for(int i = 0 ; i<ICL_N ; i++){
		//y_cl.T*(Mb - y_cl_theta)	
		/*
		ICL_control_term, 8, 1		//x
		y_cltheta ,3, ,1			//y_cl*theta
		Mb_y_cltheta ,3 ,1			//Mb-y_cl*theta
		y_clT ,8 ,3					//y_cl.transpose
		y_clT_Mb_y_cltheta ,8 ,1	//y_clT*Mb-y_cl*theta
		*/
		//y_cl_theta
		MAT_MULT(&sigma_array[i].y_cl, &theta, &y_cltheta);
		//Mb - y_cl_theta
		MAT_SUB(&sigma_array[i].Mb, &y_cltheta, &Mb_y_cltheta);
		//y_cl.T*(Mb - y_cl_theta)	
		MAT_TRANS(&sigma_array[i].y_cl,&y_clT);
		MAT_MULT(&y_clT, &Mb_y_cltheta, &y_clT_Mb_y_cltheta);

		mat_data(ICL_control_term)[0] = mat_data(ICL_control_term)[0] + mat_data(y_clT_Mb_y_cltheta)[0]; 
		mat_data(ICL_control_term)[1] = mat_data(ICL_control_term)[1] + mat_data(y_clT_Mb_y_cltheta)[1];
		mat_data(ICL_control_term)[2] = mat_data(ICL_control_term)[2] + mat_data(y_clT_Mb_y_cltheta)[2];
		mat_data(ICL_control_term)[3] = mat_data(ICL_control_term)[3] + mat_data(y_clT_Mb_y_cltheta)[3];
		mat_data(ICL_control_term)[4] = mat_data(ICL_control_term)[4] + mat_data(y_clT_Mb_y_cltheta)[4];
		mat_data(ICL_control_term)[5] = mat_data(ICL_control_term)[5] + mat_data(y_clT_Mb_y_cltheta)[5];
		mat_data(ICL_control_term)[6] = mat_data(ICL_control_term)[6] + mat_data(y_clT_Mb_y_cltheta)[6];
		mat_data(ICL_control_term)[7] = mat_data(ICL_control_term)[7] + mat_data(y_clT_Mb_y_cltheta)[7];
	}
	mat_data(ICL_theta_hat_dot)[0] = adaptive_gamma_k_icl[0]*mat_data(ICL_control_term)[0];
	mat_data(ICL_theta_hat_dot)[1] = adaptive_gamma_k_icl[1]*mat_data(ICL_control_term)[1];
	mat_data(ICL_theta_hat_dot)[2] = adaptive_gamma_k_icl[2]*mat_data(ICL_control_term)[2];
	mat_data(ICL_theta_hat_dot)[3] = adaptive_gamma_k_icl[3]*mat_data(ICL_control_term)[3];
	mat_data(ICL_theta_hat_dot)[4] = adaptive_gamma_k_icl[4]*mat_data(ICL_control_term)[4];
	mat_data(ICL_theta_hat_dot)[5] = adaptive_gamma_k_icl[5]*mat_data(ICL_control_term)[5];
	mat_data(ICL_theta_hat_dot)[6] = adaptive_gamma_k_icl[6]*mat_data(ICL_control_term)[6];
	mat_data(ICL_theta_hat_dot)[7] = adaptive_gamma_k_icl[7]*mat_data(ICL_control_term)[7];

	//theta update
	mat_data(theta)[0] = mat_data(theta)[0] + mat_data(adaptive_theta_hat_dot)[0] + mat_data(ICL_theta_hat_dot)[0];
	mat_data(theta)[1] = mat_data(theta)[1] + mat_data(adaptive_theta_hat_dot)[1] + mat_data(ICL_theta_hat_dot)[1];
	mat_data(theta)[2] = mat_data(theta)[2] + mat_data(adaptive_theta_hat_dot)[2] + mat_data(ICL_theta_hat_dot)[2];
	mat_data(theta)[3] = mat_data(theta)[3] + mat_data(adaptive_theta_hat_dot)[3] + mat_data(ICL_theta_hat_dot)[3];
	mat_data(theta)[4] = mat_data(theta)[4] + mat_data(adaptive_theta_hat_dot)[4] + mat_data(ICL_theta_hat_dot)[4];
	mat_data(theta)[5] = mat_data(theta)[5] + mat_data(adaptive_theta_hat_dot)[5] + mat_data(ICL_theta_hat_dot)[5];
	mat_data(theta)[6] = mat_data(theta)[6] + mat_data(adaptive_theta_hat_dot)[6] + mat_data(ICL_theta_hat_dot)[6];
	mat_data(theta)[7] = mat_data(theta)[7] + mat_data(adaptive_theta_hat_dot)[7] + mat_data(ICL_theta_hat_dot)[7];

	//Y_theta
	MAT_MULT(&Y, &theta, &Ytheta);
	
	/* control input M1, M2, M3 */
	output_moments[0] = -krx*mat_data(eR)[0] -kwx*mat_data(eW)[0] +mat_data(Ytheta)[0]; 
	output_moments[1] = -kry*mat_data(eR)[1] -kwy*mat_data(eW)[1] +mat_data(Ytheta)[1];
	output_moments[2] = -krz*mat_data(eR)[2] -kwz*mat_data(eW)[2] +mat_data(Ytheta)[2];
	
	mat_data(W_last5)[0] = mat_data(W_last4)[0];
	mat_data(W_last5)[1] = mat_data(W_last4)[1];
	mat_data(W_last5)[2] = mat_data(W_last4)[2];
	mat_data(W_last4)[0] = mat_data(W_last3)[0];
	mat_data(W_last4)[1] = mat_data(W_last3)[1];
	mat_data(W_last4)[2] = mat_data(W_last3)[2];
	mat_data(W_last3)[0] = mat_data(W_last2)[0];
	mat_data(W_last3)[1] = mat_data(W_last2)[1];
	mat_data(W_last3)[2] = mat_data(W_last2)[2];
	mat_data(W_last2)[0] = mat_data(W_last)[0];
	mat_data(W_last2)[1] = mat_data(W_last)[1];
	mat_data(W_last2)[2] = mat_data(W_last)[2];
	mat_data(W_last)[0] = mat_data(W)[0];
	mat_data(W_last)[1] = mat_data(W)[1];
	mat_data(W_last)[2] = mat_data(W)[2];
	mat_data(M_last5)[0] = mat_data(M_last4)[0];
	mat_data(M_last5)[1] = mat_data(M_last4)[1];
	mat_data(M_last5)[2] = mat_data(M_last4)[2];
	mat_data(M_last4)[0] = mat_data(M_last3)[0];
	mat_data(M_last4)[1] = mat_data(M_last3)[1];
	mat_data(M_last4)[2] = mat_data(M_last3)[2];
	mat_data(M_last3)[0] = mat_data(M_last2)[0];
	mat_data(M_last3)[1] = mat_data(M_last2)[1];
	mat_data(M_last3)[2] = mat_data(M_last2)[2];
	mat_data(M_last2)[0] = mat_data(M_last)[0];
	mat_data(M_last2)[1] = mat_data(M_last)[1];
	mat_data(M_last2)[2] = mat_data(M_last)[2];
	mat_data(M_last)[0] = output_moments[0];
	mat_data(M_last)[1] = output_moments[1];
	mat_data(M_last)[2] = output_moments[2];
	output_force_last5 = output_force_last4;
	output_force_last4 = output_force_last3;
	output_force_last3 = output_force_last2;
	output_force_last2 = output_force_last;
	output_force_last = *output_force;
	
#endif

}

#define l_div_4 (0.25f * (1.0f / MOTOR_TO_CG_LENGTH_M))
#define b_div_4 (+0.25f * (1.0f / COEFFICIENT_YAW))
void mr_geometry_ctrl_thrust_allocation(float *moment, float total_force)
{
	/* quadrotor thrust allocation */
	float distributed_force = total_force *= 0.25; //split force to 4 motors
	float motor_force[4];
	motor_force[0] = -l_div_4 * moment[0] + l_div_4 * moment[1] +
	                 -b_div_4 * moment[2] + distributed_force;
	motor_force[1] = +l_div_4 * moment[0] + l_div_4 * moment[1] +
	                 +b_div_4 * moment[2] + distributed_force;
	motor_force[2] = +l_div_4 * moment[0] - l_div_4 * moment[1] +
	                 -b_div_4 * moment[2] + distributed_force;
	motor_force[3] = -l_div_4 * moment[0] - l_div_4 * moment[1] +
	                 +b_div_4 * moment[2] + distributed_force;

	set_motor_value(MOTOR1, convert_motor_thrust_to_cmd(motor_force[0]));
	set_motor_value(MOTOR2, convert_motor_thrust_to_cmd(motor_force[1]));
	set_motor_value(MOTOR3, convert_motor_thrust_to_cmd(motor_force[2]));
	set_motor_value(MOTOR4, convert_motor_thrust_to_cmd(motor_force[3]));
}

void rc_mode_handler_geometry_ctrl(radio_t *rc)
{
	static bool auto_flight_mode_last = false;

	multirotor_rc_special_function_handler(rc);

	if(rc->safety == true) {
		if(rc->auto_flight == true) {
			autopilot_set_mode(AUTOPILOT_HOVERING_MODE);
		} else {
			autopilot_set_mode(AUTOPILOT_MANUAL_FLIGHT_MODE);
		}
	}

	//if mode switched to auto-flight
	if(rc->auto_flight == true && auto_flight_mode_last != true) {
		autopilot_set_mode(AUTOPILOT_HOVERING_MODE);

		//set desired position to current position
		float curr_pos[3] = {0.0f};
		get_enu_position(curr_pos);
		autopilot_assign_pos_target(curr_pos[0], curr_pos[1], curr_pos[2]);
		autopilot_assign_zero_vel_target();      //set desired velocity to zero
		autopilot_assign_zero_acc_feedforward(); //set acceleration feedforward to zero

		reset_geometry_tracking_error_integral();
	}

	if(rc->auto_flight == false) {
		autopilot_set_mode(AUTOPILOT_MANUAL_FLIGHT_MODE);
		autopilot_mission_reset();

		autopilot_assign_pos_target(0.0f, 0.0f, 0.0f);
		autopilot_assign_zero_vel_target();
		autopilot_assign_zero_acc_feedforward();

		reset_geometry_tracking_error_integral();
	}

	auto_flight_mode_last = rc->auto_flight;
}

void multirotor_geometry_control(radio_t *rc, float *desired_heading)
{
	/* check rc events */
	rc_mode_handler_geometry_ctrl(rc);

	/* get sensor status */
	bool xy_pos_available = is_xy_position_info_available();
	bool height_availabe = is_height_info_available();
	bool heading_available = is_compass_available();

	/* get imu datay */
	float accel_lpf[3];
	float gyro_lpf[3];
	get_accel_lpf(accel_lpf);
	get_gyro_lpf(gyro_lpf);

	/* get attitude quaternion */
	float attitude_q[4];
	get_attitude_quaternion(attitude_q);

	/* get roll, pitch, yaw angles */
	float attitude_roll, attitude_pitch, attitude_yaw;
	get_attitude_euler_angles(&attitude_roll, &attitude_pitch, &attitude_yaw);

	/* get direction consine matrix of current attitude */
	get_rotation_matrix_b2i(&mat_data(R));
	get_rotation_matrix_i2b(&mat_data(Rt));

	/* prepare position and velocity data */
	float curr_pos_enu[3] = {0.0f}, curr_pos_ned[3] = {0.0f};
	float curr_vel_enu[3] = {0.0f}, curr_vel_ned[3] = {0.0f};
	get_enu_position(curr_pos_enu);
	get_enu_velocity(curr_vel_enu);
	assign_vector_3x1_enu_to_ned(curr_pos_ned, curr_pos_enu);
	assign_vector_3x1_enu_to_ned(curr_vel_ned, curr_vel_enu);

	/* prepare gyroscope data */
	float gyro[3] = {0.0};
	gyro[0] = deg_to_rad(gyro_lpf[0]);
	gyro[1] = deg_to_rad(gyro_lpf[1]);
	gyro[2] = deg_to_rad(gyro_lpf[2]);

	/* prepare manual control attitude commands (euler angle) */
	euler_t attitude_cmd;
	attitude_cmd.roll = deg_to_rad(-rc->roll);
	attitude_cmd.pitch = deg_to_rad(-rc->pitch);
	if(heading_available == true) {
		//yaw control mode
		attitude_cmd.yaw = deg_to_rad(*desired_heading);
	} else {
		//yaw rate control mode
		attitude_cmd.yaw = deg_to_rad(-rc->yaw);
	}

	/* guidance loop (autopilot) */
	autopilot_guidance_handler(curr_pos_enu, curr_vel_enu);

	/* prepare desired position, velocity and acceleration feedforward */
	float pos_des_enu[3], vel_des_enu[3], accel_ff_enu[3];
	autopilot_get_pos_setpoint(pos_des_enu);
	autopilot_get_vel_setpoint(vel_des_enu);
	autopilot_get_accel_feedforward(accel_ff_enu);

	float control_moments[3] = {0.0f}, control_force = 0.0f;

	if(rc->auto_flight == true && height_availabe && heading_available) {
		if(xy_pos_available == false) {
			height_ctrl_only = true;
		}

		/* auto-flight mode (position, velocity and attitude control) */
		geometry_tracking_ctrl(&attitude_cmd, attitude_q, gyro,
		                       pos_des_enu, vel_des_enu, accel_ff_enu,
		                       curr_pos_ned, curr_vel_ned, control_moments,
		                       &control_force, height_ctrl_only);
	} else {
		/* manual flight mode (attitude control only) */
		geometry_manual_ctrl(&attitude_cmd, attitude_q, gyro, control_moments,
		                     heading_available);

		/* generate total thrust for quadrotor (open-loop) */
		control_force = 4.0f * convert_motor_cmd_to_thrust(rc->throttle * 0.01 /* [%] */);
	}

	if(rc->safety == true) {
		*desired_heading = attitude_yaw;
		barometer_set_sea_level();
		set_rgb_led_service_motor_lock_flag(true);
	} else {
		set_rgb_led_service_motor_lock_flag(false);
	}

	bool lock_motor = false;

	//lock motor if throttle values is lower than 10% during manual flight
	lock_motor |= check_motor_lock_condition(rc->throttle < 10.0f &&
	                autopilot_get_mode() == AUTOPILOT_MANUAL_FLIGHT_MODE);
	//lock motor if desired height is lower than threshold value in the takeoff mode
	lock_motor |= check_motor_lock_condition(pos_des_enu[2] < 0.10f &&
	                autopilot_get_mode() == AUTOPILOT_TAKEOFF_MODE);
	//lock motor if current position is very close to ground in the hovering mode
	lock_motor |= check_motor_lock_condition(curr_pos_enu[2] < 0.10f &&
	                autopilot_get_mode() == AUTOPILOT_HOVERING_MODE);
	//lock motor if motors are locked by autopilot
	lock_motor |= check_motor_lock_condition(autopilot_get_mode() == AUTOPILOT_MOTOR_LOCKED_MODE);
	//lock motor if radio safety botton is on
	lock_motor |= check_motor_lock_condition(rc->safety == true);

	if(lock_motor == false) {
		mr_geometry_ctrl_thrust_allocation(control_moments, control_force);
	} else {
		motor_halt();
	}
}

void send_geometry_moment_ctrl_debug(debug_msg_t *payload)
{
	float roll_error = rad_to_deg(mat_data(eR)[0]);
	float pitch_error = rad_to_deg(mat_data(eR)[1]);
	float yaw_error = rad_to_deg(mat_data(eR)[2]);

	float wx_error = rad_to_deg(mat_data(eW)[0]);
	float wy_error = rad_to_deg(mat_data(eW)[1]);
	float wz_error = rad_to_deg(mat_data(eW)[2]);

	float geometry_ctrl_feedback_moments[3];
	float geometry_ctrl_feedfoward_moments[3];

	/* calculate the feedback moment and convert the unit from [gram force * m] to [newton * m] */
	geometry_ctrl_feedback_moments[0] = (-krx*mat_data(eR)[0] -kwx*mat_data(eW)[0]);
	geometry_ctrl_feedback_moments[1] = (-kry*mat_data(eR)[1] -kwy*mat_data(eW)[1]);
	geometry_ctrl_feedback_moments[2] = (-krz*mat_data(eR)[2] -kwz*mat_data(eW)[2]);

	geometry_ctrl_feedfoward_moments[0] = mat_data(inertia_effect)[0];
	geometry_ctrl_feedfoward_moments[1] = mat_data(inertia_effect)[1];
	geometry_ctrl_feedfoward_moments[2] = mat_data(inertia_effect)[2];

	pack_debug_debug_message_header(payload, MESSAGE_ID_GEOMETRY_MOMENT_CTRL);
	pack_debug_debug_message_float(&roll_error, payload);
	pack_debug_debug_message_float(&pitch_error, payload);
	pack_debug_debug_message_float(&yaw_error, payload);
	pack_debug_debug_message_float(&wx_error, payload);
	pack_debug_debug_message_float(&wy_error, payload);
	pack_debug_debug_message_float(&wz_error, payload);
	pack_debug_debug_message_float(&geometry_ctrl_feedback_moments[0], payload);
	pack_debug_debug_message_float(&geometry_ctrl_feedback_moments[1], payload);
	pack_debug_debug_message_float(&geometry_ctrl_feedback_moments[2], payload);
	pack_debug_debug_message_float(&geometry_ctrl_feedfoward_moments[0], payload);
	pack_debug_debug_message_float(&geometry_ctrl_feedfoward_moments[1], payload);
	pack_debug_debug_message_float(&geometry_ctrl_feedfoward_moments[2], payload);
}

void send_controller_estimation_adaptive_debug(debug_msg_t *payload)
{
	//y_cl_theta
	MAT_ALLOC(y_cltheta_test,3,1);
	MAT_INIT(y_cltheta_test,3,1);
	MAT_MULT(&sigma_array[ICL_sigma_index].y_cl, &theta, &y_cltheta_test);
	//Mb - y_cl_theta
	MAT_ALLOC(Mb_y_cltheta_test,3,1);
	MAT_INIT(Mb_y_cltheta_test,3,1);
	MAT_SUB(&sigma_array[ICL_sigma_index].Mb, &y_cltheta_test, &Mb_y_cltheta_test);
	//y_cl.T*(Mb - y_cl_theta)	
	float roll_error = rad_to_deg(mat_data(eR)[0]);
	float pitch_error = rad_to_deg(mat_data(eR)[1]);
	float yaw_error = rad_to_deg(mat_data(eR)[2]);
	
	float theta0 = mat_data(theta)[0];
	float theta1 = mat_data(theta)[1];
	float theta2 = mat_data(theta)[2];
	float theta3 = mat_data(theta)[3];
	float theta4 = mat_data(theta)[4];
	float theta5 = mat_data(theta)[5];
	float theta6 = mat_data(theta)[6];
	float theta7 = mat_data(theta)[7];
	
	pack_debug_debug_message_header(payload, MESSAGE_ID_ADAPTIVE_THETA);
	pack_debug_debug_message_float(&roll_error, payload);
	pack_debug_debug_message_float(&pitch_error, payload);
	pack_debug_debug_message_float(&yaw_error, payload);
	
	pack_debug_debug_message_float(&theta0, payload);
	pack_debug_debug_message_float(&theta1, payload);
	pack_debug_debug_message_float(&theta2, payload);
	pack_debug_debug_message_float(&theta3, payload);
	pack_debug_debug_message_float(&theta4, payload);
	pack_debug_debug_message_float(&theta5, payload);
	pack_debug_debug_message_float(&theta6, payload);
	pack_debug_debug_message_float(&theta7, payload);
}
void send_geometry_tracking_ctrl_debug(debug_msg_t *payload)
{
	pack_debug_debug_message_header(payload, MESSAGE_ID_GEOMETRY_TRACKING_CTRL);
	pack_debug_debug_message_float(&pos_error[0], payload);
	pack_debug_debug_message_float(&pos_error[1], payload);
	pack_debug_debug_message_float(&pos_error[2], payload);
}

void send_uav_dynamics_debug(debug_msg_t *payload)
{
	pack_debug_debug_message_header(payload, MESSAGE_ID_UAV_DYNAMICS_DEBUG);
	pack_debug_debug_message_float(&uav_dynamics_m[0], payload);
	pack_debug_debug_message_float(&uav_dynamics_m[1], payload);
	pack_debug_debug_message_float(&uav_dynamics_m[2], payload);
	pack_debug_debug_message_float(&uav_dynamics_m_rot_frame[0], payload);
	pack_debug_debug_message_float(&uav_dynamics_m_rot_frame[1], payload);
	pack_debug_debug_message_float(&uav_dynamics_m_rot_frame[2], payload);
}
