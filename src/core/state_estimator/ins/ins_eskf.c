#include <math.h>
#include "arm_math.h"
#include "ahrs.h"
#include "lpf.h"
#include "matrix.h"
#include "se3_math.h"
#include "quaternion.h"

/*==================================================================*
 * error-state kalman filter for multirotor full state estimation   *
 * derived by Sheng-Wen, Cheng (shengwen1997.tw@gmail.com)          *
 * estimate state = [position(3x1); velocity(3x1); quaternion(4x1)] *
 *==================================================================*/

#define ESKF_RESCALE(number) (number * 10e7) //to improve the numerical stability

#define R(r, c)             _R.pData[(r * 3) + c]
#define Rt(r, c)            _Rt.pData[(r * 3) + c]
#define P_prior(r, c)       _P_prior.pData[(r * 9) + c]
#define P_post(r, c)        _P_post.pData[(r * 9) + c]
#define R_am_ab_dt(r, c)    _R_am_ab_dt.pData[(r * 3) + c]
#define Rt_wm_wb_dt(r, c)   _Rt_wm_wb_dt.pData[(r * 3) + c]
#define Q_i(r, c)           _Q_i.pData[(r * 6) + c]
#define K_accel(r, c)       _K_accel.pData[(r * 3) + c]
#define K_mag(r, c)         _K_mag.pData[(r * 3) + c]
#define K_gps(r, c)         _K_gps.pData[(r * 4) + c]
#define K_baro(r, c)        _K_baro.pData[(r * 2) + c]
#define V_accel(r, c)       _V_accel.pData[(r * 3) + c]
#define V_mag(r, c)         _V_mag.pData[(r * 3) + c]
#define PHt_accel(r, c)     _PHt_accel.pData[(r * 3) + c]
#define PHt_mag(r, c)       _PHt_mag.pData[(r * 3) + c]
#define PHt_gps             _PHt_gps.pData[(r * 4) + c]
#define PHt_baro            _PHt_baro.pData[(r * 2) + c]
#define HPHt_V_accel(r, c)  _HPHt_V_accel.pData[(r * 3) + c]
#define HPHt_V_mag(r, c)    _HPHt_V_mag.pData[(r * 3) + c]
#define HPHt_V_gps          _HPHt_V_gps.pData[(r * 4) + c]
#define HPHt_V_baro         _HPHt_V_baro.pData[(r * 2) + c]

MAT_ALLOC(x_nominal, 10, 1);
MAT_ALLOC(x_error_state, 9, 1);
MAT_ALLOC(_Q_i, 6, 6);
MAT_ALLOC(_V_accel, 3, 3);
MAT_ALLOC(_V_mag, 3, 3);
MAT_ALLOC(_P_prior, 9, 9);
MAT_ALLOC(_P_post, 9, 9);
MAT_ALLOC(_R, 3, 3);
MAT_ALLOC(_Rt, 3, 3);
MAT_ALLOC(_R_am_ab_dt, 3, 3);
MAT_ALLOC(_Rt_wm_wb_dt, 3, 3);
MAT_ALLOC(_PHt_accel, 9, 3);
MAT_ALLOC(_PHt_mag, 9, 3);
MAT_ALLOC(_PHt_gps, 9, 4);
MAT_ALLOC(_PHt_baro, 9, 2);
MAT_ALLOC(_HPHt_V_accel, 3, 3);
MAT_ALLOC(_HPHt_V_mag, 3, 3);
MAT_ALLOC(_HPHt_V_gps, 4, 4);
MAT_ALLOC(_HPHt_V_baro, 2, 2);
MAT_ALLOC(_HPHt_V_accel_inv, 3, 3);
MAT_ALLOC(_HPHt_V_mag_inv, 3, 3);
MAT_ALLOC(_HPHt_V_gps_inv, 4, 4);
MAT_ALLOC(_HPHt_V_baro_inv, 2, 2);
MAT_ALLOC(_K_accel, 9, 3);
MAT_ALLOC(_K_mag, 9, 3);
MAT_ALLOC(_K_gps, 9, 4);
MAT_ALLOC(_K_baro, 9, 2);

float dt;
float neg_half_dt;
float half_dt_squared;

void eskf_ins_init(float dt)
{
	dt = dt;
	neg_half_dt = -0.5f * dt;
	half_dt_squared = 0.5f * dt * dt;

	//TODO: initialize all matrices!

	/* initialize the nominal state */
	mat_data(x_nominal)[0] = 0.0f; //px
	mat_data(x_nominal)[1] = 0.0f; //py
	mat_data(x_nominal)[2] = 0.0f; //pz
	mat_data(x_nominal)[3] = 0.0f; //vx
	mat_data(x_nominal)[4] = 0.0f; //vy
	mat_data(x_nominal)[5] = 0.0f; //vz
	mat_data(x_nominal)[6] = 1.0f; //q0
	mat_data(x_nominal)[7] = 0.0f; //q1
	mat_data(x_nominal)[8] = 0.0f; //q2
	mat_data(x_nominal)[9] = 0.0f; //q3

	/* initialize _Q_i matrix */
	matrix_reset(mat_data(_Q_i), 6, 6);
	Q_i(0, 0) = ESKF_RESCALE(1e-5); //Var(ax)
	Q_i(1, 1) = ESKF_RESCALE(1e-5); //Var(ay)
	Q_i(2, 2) = ESKF_RESCALE(1e-5); //Var(az)
	Q_i(3, 3) = ESKF_RESCALE(1e-5); //Var(wx)
	Q_i(4, 4) = ESKF_RESCALE(1e-5); //Var(wy)
	Q_i(5, 5) = ESKF_RESCALE(1e-5); //Var(wz)

	/* initialize P matrix */
	matrix_reset(mat_data(_P_post), 9, 9);
	P_post(0, 0) = ESKF_RESCALE(5.0f); //Var(px)
	P_post(1, 1) = ESKF_RESCALE(5.0f); //Var(py)
	P_post(2, 2) = ESKF_RESCALE(5.0f); //Var(pz)
	P_post(3, 3) = ESKF_RESCALE(5.0f); //Var(vx)
	P_post(4, 4) = ESKF_RESCALE(5.0f); //Var(vy)
	P_post(5, 5) = ESKF_RESCALE(5.0f); //Var(vz)
	P_post(6, 6) = ESKF_RESCALE(5.0f); //Var(theta_x)
	P_post(7, 7) = ESKF_RESCALE(5.0f); //Var(theta_y)
	P_post(8, 8) = ESKF_RESCALE(5.0f); //Var(theta_z)

	/* initialize V_accel matrix */
	matrix_reset(mat_data(_V_accel), 9, 9);
	V_accel(0, 0) = ESKF_RESCALE(7e-1); //Var(gx)
	V_accel(1, 1) = ESKF_RESCALE(7e-1); //Var(gy)
	V_accel(2, 2) = ESKF_RESCALE(7e-1); //Var(gz)

	/* initialize V_mag matrix */
	matrix_reset(mat_data(_V_mag), 9, 9);
	V_mag(0, 0) = ESKF_RESCALE(3); //Var(mx)
	V_mag(1, 1) = ESKF_RESCALE(3); //Var(my)
	V_mag(2, 2) = ESKF_RESCALE(3); //Var(mz)

	/* initial V_gps matrix */
	//TODO

	/* initialize V_barometer matrix */
	//TODO
}

void eskf_ins_predict(float *accel, float *gyro)
{
	/* input variables (ned frame) */
	float accel_b_x = accel[0];
	float accel_b_y = accel[1];
	float accel_b_z = accel[2];
	float gyro_b_x = gyro[0];
	float gyro_b_y = gyro[1];
	float gyro_b_z = gyro[2];

	/* body-frame to inertial-frame conversion */
	float accel_i_ned[3] = {0};
	accel_i_ned[0] = Rt(0, 0) * accel_b_x + Rt(0, 1) * accel_b_y + Rt(0, 2) * accel_b_z;
	accel_i_ned[1] = Rt(1, 0) * accel_b_x + Rt(1, 1) * accel_b_y + Rt(1, 2) * accel_b_z;
	accel_i_ned[2] = Rt(2, 0) * accel_b_x + Rt(2, 1) * accel_b_y + Rt(1, 2) * accel_b_z;

	//ned to enu conversion
	float accel_i[3];
	accel_i[0] =  accel_i_ned[1];
	accel_i[1] =  accel_i_ned[0];
	accel_i[2] = -accel_i_ned[2];

	/*======================*
	 * nominal state update *
	 *======================*/

	//velocity integration
	mat_data(x_nominal)[0] += accel_i[0] * dt;
	mat_data(x_nominal)[1] += accel_i[1] * dt;
	mat_data(x_nominal)[2] += accel_i[2] * dt;

	//position integration
	mat_data(x_nominal)[0] += (mat_data(x_nominal)[3] * dt) +
	                          (accel_i[0] * half_dt_squared);
	mat_data(x_nominal)[1] += (mat_data(x_nominal)[4] * dt) +
	                          (accel_i[1] * half_dt_squared);
	mat_data(x_nominal)[2] += (mat_data(x_nominal)[5] * dt) +
	                          (accel_i[2] * half_dt_squared);

	/* calculate quaternion time derivative */
	float w[4];
	w[0] = 0.0f;
	w[1] = gyro[0];
	w[2] = gyro[1];
	w[3] = gyro[2];
	float q_dot[4];
	quaternion_mult(w, &mat_data(x_nominal)[6], q_dot);

	//quaternion integration
	mat_data(x_nominal)[6] = mat_data(x_nominal)[6] + (q_dot[0] * neg_half_dt);
	mat_data(x_nominal)[7] = mat_data(x_nominal)[7] + (q_dot[1] * neg_half_dt);
	mat_data(x_nominal)[8] = mat_data(x_nominal)[8] + (q_dot[2] * neg_half_dt);
	mat_data(x_nominal)[9] = mat_data(x_nominal)[9] + (q_dot[3] * neg_half_dt);
	quat_normalize(mat_data(x_nominal));

	/*==================================*
	 * process covatiance matrix update *
	 *==================================*/

	/* codeblock for preventing nameing conflict */
	{
		/* -R * skew_matrix(a_m - a_b) * dt */
		R_am_ab_dt(0, 0) = -dt*(R(0,1)*accel_b_z-R(0,2)*accel_b_y);
		R_am_ab_dt(0, 1) = dt*(R(0,0)*accel_b_z-R(0,2)*accel_b_x);
		R_am_ab_dt(0, 2) = -dt*(R(0,0)*accel_b_y-R(0,1)*accel_b_x);
		R_am_ab_dt(1, 0) = -dt*(R(1,1)*accel_b_z-R(1,2)*accel_b_y);
		R_am_ab_dt(1, 1) = dt*(R(1,0)*accel_b_z-R(1,2)*accel_b_x);
		R_am_ab_dt(1, 2) = -dt*(R(1,0)*accel_b_y-R(1,1)*accel_b_x);
		R_am_ab_dt(2, 0) = -dt*(R(2,1)*accel_b_z-R(2,2)*accel_b_y);
		R_am_ab_dt(2, 1) = dt*(R(2,0)*accel_b_z-R(2,2)*accel_b_x);
		R_am_ab_dt(2, 2) = -dt*(R(2,0)*accel_b_y-R(2,1)*accel_b_x);

		/* Rt{w_m - w_b} * dt */
		float c0 = dt*gyro_b_z;
		float c1 = dt*gyro_b_y;
		float c2 = dt*gyro_b_x;

		Rt_wm_wb_dt(0, 0) = 1.0;
		Rt_wm_wb_dt(0, 1) = c0;
		Rt_wm_wb_dt(0, 2) = -c1;
		Rt_wm_wb_dt(1, 0) = -c0;
		Rt_wm_wb_dt(1, 1) = 1.0;
		Rt_wm_wb_dt(1, 2) = c2;
		Rt_wm_wb_dt(2, 0) = c1;
		Rt_wm_wb_dt(2, 1) = -c2;
		Rt_wm_wb_dt(2, 2) = 1.0;
	}

	/* codeblock for preventing nameing conflict */
	{
		/* calculate a priori process covariance matrix */
		float c0 = P_post(5,8)+P_post(6,8)*R_am_ab_dt(2,0)+P_post(7,8)*R_am_ab_dt(2,1)+P_post(8,8)*R_am_ab_dt(2,2);
		float c1 = P_post(5,7)+P_post(6,7)*R_am_ab_dt(2,0)+P_post(7,7)*R_am_ab_dt(2,1)+P_post(8,7)*R_am_ab_dt(2,2);
		float c2 = P_post(5,6)+P_post(6,6)*R_am_ab_dt(2,0)+P_post(7,6)*R_am_ab_dt(2,1)+P_post(8,6)*R_am_ab_dt(2,2);
		float c3 = P_post(4,8)+P_post(6,8)*R_am_ab_dt(1,0)+P_post(7,8)*R_am_ab_dt(1,1)+P_post(8,8)*R_am_ab_dt(1,2);
		float c4 = P_post(4,7)+P_post(6,7)*R_am_ab_dt(1,0)+P_post(7,7)*R_am_ab_dt(1,1)+P_post(8,7)*R_am_ab_dt(1,2);
		float c5 = P_post(4,6)+P_post(6,6)*R_am_ab_dt(1,0)+P_post(7,6)*R_am_ab_dt(1,1)+P_post(8,6)*R_am_ab_dt(1,2);
		float c6 = P_post(3,8)+P_post(6,8)*R_am_ab_dt(0,0)+P_post(7,8)*R_am_ab_dt(0,1)+P_post(8,8)*R_am_ab_dt(0,2);
		float c7 = P_post(3,7)+P_post(6,7)*R_am_ab_dt(0,0)+P_post(7,7)*R_am_ab_dt(0,1)+P_post(8,7)*R_am_ab_dt(0,2);
		float c8 = P_post(3,6)+P_post(6,6)*R_am_ab_dt(0,0)+P_post(7,6)*R_am_ab_dt(0,1)+P_post(8,6)*R_am_ab_dt(0,2);
		float c9 = P_post(6,8)*Rt_wm_wb_dt(2,0)+P_post(7,8)*Rt_wm_wb_dt(2,1)+P_post(8,8)*Rt_wm_wb_dt(2,2);
		float c10 = P_post(6,7)*Rt_wm_wb_dt(2,0)+P_post(7,7)*Rt_wm_wb_dt(2,1)+P_post(8,7)*Rt_wm_wb_dt(2,2);
		float c11 = P_post(6,6)*Rt_wm_wb_dt(2,0)+P_post(7,6)*Rt_wm_wb_dt(2,1)+P_post(8,6)*Rt_wm_wb_dt(2,2);
		float c12 = P_post(6,8)*Rt_wm_wb_dt(1,0)+P_post(7,8)*Rt_wm_wb_dt(1,1)+P_post(8,8)*Rt_wm_wb_dt(1,2);
		float c13 = P_post(6,7)*Rt_wm_wb_dt(1,0)+P_post(7,7)*Rt_wm_wb_dt(1,1)+P_post(8,7)*Rt_wm_wb_dt(1,2);
		float c14 = P_post(6,6)*Rt_wm_wb_dt(1,0)+P_post(7,6)*Rt_wm_wb_dt(1,1)+P_post(8,6)*Rt_wm_wb_dt(1,2);
		float c15 = P_post(6,8)*Rt_wm_wb_dt(0,0)+P_post(7,8)*Rt_wm_wb_dt(0,1)+P_post(8,8)*Rt_wm_wb_dt(0,2);
		float c16 = P_post(6,7)*Rt_wm_wb_dt(0,0)+P_post(7,7)*Rt_wm_wb_dt(0,1)+P_post(8,7)*Rt_wm_wb_dt(0,2);
		float c17 = P_post(6,6)*Rt_wm_wb_dt(0,0)+P_post(7,6)*Rt_wm_wb_dt(0,1)+P_post(8,6)*Rt_wm_wb_dt(0,2);
		float c18 = P_post(2,8)+P_post(5,8)*dt;
		float c19 = P_post(2,7)+P_post(5,7)*dt;
		float c20 = P_post(2,6)+P_post(5,6)*dt;
		float c21 = P_post(1,8)+P_post(4,8)*dt;
		float c22 = P_post(1,7)+P_post(4,7)*dt;
		float c23 = P_post(1,6)+P_post(4,6)*dt;
		float c24 = P_post(0,8)+P_post(3,8)*dt;
		float c25 = P_post(0,7)+P_post(3,7)*dt;
		float c26 = P_post(0,6)+P_post(3,6)*dt;
		float c27 = P_post(5,5)*dt;
		float c28 = P_post(5,4)*dt;
		float c29 = P_post(5,3)*dt;
		float c30 = P_post(4,5)*dt;
		float c31 = P_post(4,4)*dt;
		float c32 = P_post(4,3)*dt;
		float c33 = P_post(3,5)*dt;
		float c34 = P_post(3,4)*dt;
		float c35 = P_post(3,3)*dt;
		float c63 = P_post(8,5)*R_am_ab_dt(2,2);
		float c66 = P_post(8,5)*R_am_ab_dt(1,2);
		float c67 = P_post(8,4)*R_am_ab_dt(1,2);
		float c68 = P_post(7,5)*R_am_ab_dt(2,1);
		float c72 = P_post(8,5)*R_am_ab_dt(0,2);
		float c73 = P_post(8,4)*R_am_ab_dt(0,2);
		float c74 = P_post(7,5)*R_am_ab_dt(1,1);
		float c75 = P_post(8,3)*R_am_ab_dt(0,2);
		float c76 = P_post(7,4)*R_am_ab_dt(1,1);
		float c77 = P_post(6,5)*R_am_ab_dt(2,0);
		float c81 = P_post(7,5)*R_am_ab_dt(0,1);
		float c82 = P_post(7,4)*R_am_ab_dt(0,1);
		float c83 = P_post(6,5)*R_am_ab_dt(1,0);
		float c84 = P_post(7,3)*R_am_ab_dt(0,1);
		float c85 = P_post(6,4)*R_am_ab_dt(1,0);
		float c87 = P_post(6,5)*R_am_ab_dt(0,0);
		float c88 = P_post(6,4)*R_am_ab_dt(0,0);
		float c89 = P_post(6,3)*R_am_ab_dt(0,0);
		float c90 = P_post(3,3)+c75+c84+c89;
		float c91 = P_post(3,4)+c73+c82+c88;
		float c93 = P_post(3,5)+c72+c81+c87;
		float c94 = P_post(4,4)+c67+c76+c85;
		float c96 = P_post(4,5)+c66+c74+c83;
		float c98 = P_post(5,5)+c63+c68+c77;
		float c108 = P_post(2,5)+c27;
		float c109 = P_post(2,4)+c28;
		float c110 = P_post(2,3)+c29;
		float c111 = P_post(1,5)+c30;
		float c112 = P_post(1,4)+c31;
		float c113 = P_post(1,3)+c32;
		float c114 = P_post(0,5)+c33;
		float c115 = P_post(0,4)+c34;
		float c116 = P_post(0,3)+c35;

		P_prior(0, 0) = P_post(0,0)+P_post(3,0)*dt+c116*dt;
		P_prior(0, 1) = P_post(0,1)+P_post(3,1)*dt+c115*dt;
		P_prior(0, 2) = P_post(0,2)+P_post(3,2)*dt+c114*dt;
		P_prior(0, 3) = c116+R_am_ab_dt(0,0)*c26+R_am_ab_dt(0,1)*c25+R_am_ab_dt(0,2)*c24;
		P_prior(0, 4) = c115+R_am_ab_dt(1,0)*c26+R_am_ab_dt(1,1)*c25+R_am_ab_dt(1,2)*c24;
		P_prior(0, 5) = c114+R_am_ab_dt(2,0)*c26+R_am_ab_dt(2,1)*c25+R_am_ab_dt(2,2)*c24;
		P_prior(0, 6) = Rt_wm_wb_dt(0,0)*c26+Rt_wm_wb_dt(0,1)*c25+Rt_wm_wb_dt(0,2)*c24;
		P_prior(0, 7) = Rt_wm_wb_dt(1,0)*c26+Rt_wm_wb_dt(1,1)*c25+Rt_wm_wb_dt(1,2)*c24;
		P_prior(0, 8) = Rt_wm_wb_dt(2,0)*c26+Rt_wm_wb_dt(2,1)*c25+Rt_wm_wb_dt(2,2)*c24;
		P_prior(1, 0) = P_prior(0, 1);
		P_prior(1, 1) = P_post(1,1)+P_post(4,1)*dt+c112*dt;
		P_prior(1, 2) = P_post(1,2)+P_post(4,2)*dt+c111*dt;
		P_prior(1, 3) = c113+R_am_ab_dt(0,0)*c23+R_am_ab_dt(0,1)*c22+R_am_ab_dt(0,2)*c21;
		P_prior(1, 4) = c112+R_am_ab_dt(1,0)*c23+R_am_ab_dt(1,1)*c22+R_am_ab_dt(1,2)*c21;
		P_prior(1, 5) = c111+R_am_ab_dt(2,0)*c23+R_am_ab_dt(2,1)*c22+R_am_ab_dt(2,2)*c21;
		P_prior(1, 6) = Rt_wm_wb_dt(0,0)*c23+Rt_wm_wb_dt(0,1)*c22+Rt_wm_wb_dt(0,2)*c21;
		P_prior(1, 7) = Rt_wm_wb_dt(1,0)*c23+Rt_wm_wb_dt(1,1)*c22+Rt_wm_wb_dt(1,2)*c21;
		P_prior(1, 8) = Rt_wm_wb_dt(2,0)*c23+Rt_wm_wb_dt(2,1)*c22+Rt_wm_wb_dt(2,2)*c21;
		P_prior(2, 0) = P_prior(0, 2);
		P_prior(2, 1) = P_prior(1, 2);
		P_prior(2, 2) = P_post(2,2)+P_post(5,2)*dt+c108*dt;
		P_prior(2, 3) = c110+R_am_ab_dt(0,0)*c20+R_am_ab_dt(0,1)*c19+R_am_ab_dt(0,2)*c18;
		P_prior(2, 4) = c109+R_am_ab_dt(1,0)*c20+R_am_ab_dt(1,1)*c19+R_am_ab_dt(1,2)*c18;
		P_prior(2, 5) = c108+R_am_ab_dt(2,0)*c20+R_am_ab_dt(2,1)*c19+R_am_ab_dt(2,2)*c18;
		P_prior(2, 6) = Rt_wm_wb_dt(0,0)*c20+Rt_wm_wb_dt(0,1)*c19+Rt_wm_wb_dt(0,2)*c18;
		P_prior(2, 7) = Rt_wm_wb_dt(1,0)*c20+Rt_wm_wb_dt(1,1)*c19+Rt_wm_wb_dt(1,2)*c18;
		P_prior(2, 8) = Rt_wm_wb_dt(2,0)*c20+Rt_wm_wb_dt(2,1)*c19+Rt_wm_wb_dt(2,2)*c18;
		P_prior(3, 0) = P_prior(0, 3);
		P_prior(3, 1) = P_prior(1, 3);
		P_prior(3, 2) = P_prior(2, 3);
		P_prior(3, 3) = Q_i(0,0)+c90+R_am_ab_dt(0,0)*c8+R_am_ab_dt(0,1)*c7+R_am_ab_dt(0,2)*c6;
		P_prior(3, 4) = c91+R_am_ab_dt(1,0)*c8+R_am_ab_dt(1,1)*c7+R_am_ab_dt(1,2)*c6;
		P_prior(3, 5) = c93+R_am_ab_dt(2,0)*c8+R_am_ab_dt(2,1)*c7+R_am_ab_dt(2,2)*c6;
		P_prior(3, 6) = Rt_wm_wb_dt(0,0)*c8+Rt_wm_wb_dt(0,1)*c7+Rt_wm_wb_dt(0,2)*c6;
		P_prior(3, 7) = Rt_wm_wb_dt(1,0)*c8+Rt_wm_wb_dt(1,1)*c7+Rt_wm_wb_dt(1,2)*c6;
		P_prior(3, 8) = Rt_wm_wb_dt(2,0)*c8+Rt_wm_wb_dt(2,1)*c7+Rt_wm_wb_dt(2,2)*c6;
		P_prior(4, 0) = P_prior(0, 4);
		P_prior(4, 1) = P_prior(1, 4);
		P_prior(4, 2) = P_prior(2, 4);
		P_prior(4, 3) = P_prior(3, 4);
		P_prior(4, 4) = Q_i(1,1)+c94+R_am_ab_dt(1,0)*c5+R_am_ab_dt(1,1)*c4+R_am_ab_dt(1,2)*c3;
		P_prior(4, 5) = c96+R_am_ab_dt(2,0)*c5+R_am_ab_dt(2,1)*c4+R_am_ab_dt(2,2)*c3;
		P_prior(4, 6) = Rt_wm_wb_dt(0,0)*c5+Rt_wm_wb_dt(0,1)*c4+Rt_wm_wb_dt(0,2)*c3;
		P_prior(4, 7) = Rt_wm_wb_dt(1,0)*c5+Rt_wm_wb_dt(1,1)*c4+Rt_wm_wb_dt(1,2)*c3;
		P_prior(4, 8) = Rt_wm_wb_dt(2,0)*c5+Rt_wm_wb_dt(2,1)*c4+Rt_wm_wb_dt(2,2)*c3;
		P_prior(5, 0) = P_prior(0, 5);
		P_prior(5, 1) = P_prior(1, 5);
		P_prior(5, 2) = P_prior(2, 5);
		P_prior(5, 3) = P_prior(3, 5);
		P_prior(5, 4) = P_prior(4, 5);
		P_prior(5, 5) = Q_i(2,2)+c98+R_am_ab_dt(2,0)*c2+R_am_ab_dt(2,1)*c1+R_am_ab_dt(2,2)*c0;
		P_prior(5, 6) = Rt_wm_wb_dt(0,0)*c2+Rt_wm_wb_dt(0,1)*c1+Rt_wm_wb_dt(0,2)*c0;
		P_prior(5, 7) = Rt_wm_wb_dt(1,0)*c2+Rt_wm_wb_dt(1,1)*c1+Rt_wm_wb_dt(1,2)*c0;
		P_prior(5, 8) = Rt_wm_wb_dt(2,0)*c2+Rt_wm_wb_dt(2,1)*c1+Rt_wm_wb_dt(2,2)*c0;
		P_prior(6, 0) = P_prior(0, 6);
		P_prior(6, 1) = P_prior(1, 6);
		P_prior(6, 2) = P_prior(2, 6);
		P_prior(6, 3) = P_prior(3, 6);
		P_prior(6, 4) = P_prior(4, 6);
		P_prior(6, 5) = P_prior(5, 6);
		P_prior(6, 6) = Q_i(3,3)+Rt_wm_wb_dt(0,0)*c17+Rt_wm_wb_dt(0,1)*c16+Rt_wm_wb_dt(0,2)*c15;
		P_prior(6, 7) = Rt_wm_wb_dt(1,0)*c17+Rt_wm_wb_dt(1,1)*c16+Rt_wm_wb_dt(1,2)*c15;
		P_prior(6, 8) = Rt_wm_wb_dt(2,0)*c17+Rt_wm_wb_dt(2,1)*c16+Rt_wm_wb_dt(2,2)*c15;
		P_prior(7, 0) = P_prior(0, 7);
		P_prior(7, 1) = P_prior(1, 7);
		P_prior(7, 2) = P_prior(2, 7);
		P_prior(7, 3) = P_prior(3, 7);
		P_prior(7, 4) = P_prior(4, 7);
		P_prior(7, 5) = P_prior(5, 7);
		P_prior(7, 6) = P_prior(6, 7);
		P_prior(7, 7) = Q_i(4,4)+Rt_wm_wb_dt(1,0)*c14+Rt_wm_wb_dt(1,1)*c13+Rt_wm_wb_dt(1,2)*c12;
		P_prior(7, 8) = Rt_wm_wb_dt(2,0)*c14+Rt_wm_wb_dt(2,1)*c13+Rt_wm_wb_dt(2,2)*c12;
		P_prior(8, 0) = P_prior(0, 8);
		P_prior(8, 1) = P_prior(1, 8);
		P_prior(8, 2) = P_prior(2, 8);
		P_prior(8, 3) = P_prior(3, 8);
		P_prior(8, 4) = P_prior(4, 8);
		P_prior(8, 5) = P_prior(5, 8);
		P_prior(8, 6) = P_prior(6, 8);
		P_prior(8, 7) = P_prior(7, 8);
		P_prior(8, 8) = Q_i(5,5)+Rt_wm_wb_dt(2,0)*c11+Rt_wm_wb_dt(2,1)*c10+Rt_wm_wb_dt(2,2)*c9;
	}

	/*=================================================*
	 * convert estimated quaternion to R and Rt matrix *
	 *=================================================*/
	float *q = &mat_data(x_nominal)[0];
	quat_to_rotation_matrix(q, mat_data(_Rt), mat_data(_R));
}

void eskf_ins_accelerometer_correct(float *accel)
{
	float gx = accel[0];
	float gy = accel[1];
	float gz = accel[2];

	float q0 = mat_data(x_nominal)[6];
	float q1 = mat_data(x_nominal)[7];
	float q2 = mat_data(x_nominal)[8];
	float q3 = mat_data(x_nominal)[9];

	/* codeblock for preventing nameing conflict */
	{
	}

	/* codeblock for preventing nameing conflict */
	{
	}

	/* calculate kalman gain */
	//K = P * Ht * inv(H*P*Ht + V)
	MAT_INV(&_HPHt_V_accel, &_HPHt_V_accel_inv);
	MAT_MULT(&_PHt_accel, &_HPHt_V_accel_inv, &_K_accel);

	/* codeblock for preventing nameing conflict */
	{
		/* calculate error state residual */
		float c0 = gz-q0*q0+q1*q1+q2*q2-q3*q3;
		float c1 = -gx+q0*q2*2.0+q1*q3*2.0;
		float c2 = gy+q0*q1*2.0-q2*q3*2.0;

		mat_data(x_error_state)[6] = -K_accel(6,0)*c1+K_accel(6,2)*c0+K_accel(6,1)*c2;
		mat_data(x_error_state)[7] = -K_accel(7,0)*c1+K_accel(7,2)*c0+K_accel(7,1)*c2;
		mat_data(x_error_state)[8] = -K_accel(8,0)*c1+K_accel(8,2)*c0+K_accel(8,1)*c2;
	}

	/* codeblock for preventing nameing conflict */
	{
		/* calculate a posteriori process covariance matrix */
		//P = (I - K*H) * P
		float c0_ = -q3*q3;
		float c1_ = q0*q3*2.0;
		float c2_ = q2*q2;
		float c3_ = q1*q1;
		float c4_ = q0*q0;
		float c5_ = q1*q2*2.0;

		float c0 = c0_-c2_+c3_+c4_;
		float c1 = c0_+c2_-c3_+c4_;
		float c2 = c1_-c5_;
		float c3 = q0*q2*2.0-q1*q3*2.0;
		float c4 = c1_+c5_;
		float c5 = q0*q1*2.0+q2*q3*2.0;
		float c6 = K_accel(7,0)*c0+K_accel(7,1)*c4-K_accel(7,2)*c3-1.0;
		float c7 = -K_accel(6,0)*c2+K_accel(6,1)*c1+K_accel(6,2)*c5+1.0;
		float c8 = -K_accel(8,0)*c2+K_accel(8,1)*c1+K_accel(8,2)*c5;
		float c9 = K_accel(8,0)*c0+K_accel(8,1)*c4-K_accel(8,2)*c3;
		float c10 = -K_accel(7,0)*c2+K_accel(7,1)*c1+K_accel(7,2)*c5;
		float c11 = K_accel(6,0)*c0+K_accel(6,1)*c4-K_accel(6,2)*c3;
		float c12 = -K_accel(5,0)*c2+K_accel(5,1)*c1+K_accel(5,2)*c5;
		float c13 = K_accel(5,0)*c0+K_accel(5,1)*c4-K_accel(5,2)*c3;
		float c14 = -K_accel(4,0)*c2+K_accel(4,1)*c1+K_accel(4,2)*c5;
		float c15 = K_accel(4,0)*c0+K_accel(4,1)*c4-K_accel(4,2)*c3;
		float c16 = -K_accel(3,0)*c2+K_accel(3,1)*c1+K_accel(3,2)*c5;
		float c17 = K_accel(3,0)*c0+K_accel(3,1)*c4-K_accel(3,2)*c3;
		float c18 = -K_accel(2,0)*c2+K_accel(2,1)*c1+K_accel(2,2)*c5;
		float c19 = K_accel(2,0)*c0+K_accel(2,1)*c4-K_accel(2,2)*c3;
		float c20 = -K_accel(1,0)*c2+K_accel(1,1)*c1+K_accel(1,2)*c5;
		float c21 = K_accel(1,0)*c0+K_accel(1,1)*c4-K_accel(1,2)*c3;
		float c22 = -K_accel(0,0)*c2+K_accel(0,1)*c1+K_accel(0,2)*c5;
		float c23 = K_accel(0,0)*c0+K_accel(0,1)*c4-K_accel(0,2)*c3;

		P_post(0, 0) = P_prior(0,0)+P_prior(6,0)*c22-P_prior(7,0)*c23;
		P_post(0, 1) = P_prior(0,1)+P_prior(6,1)*c22-P_prior(7,1)*c23;
		P_post(0, 2) = P_prior(0,2)+P_prior(6,2)*c22-P_prior(7,2)*c23;
		P_post(0, 3) = P_prior(0,3)+P_prior(6,3)*c22-P_prior(7,3)*c23;
		P_post(0, 4) = P_prior(0,4)+P_prior(6,4)*c22-P_prior(7,4)*c23;
		P_post(0, 5) = P_prior(0,5)+P_prior(6,5)*c22-P_prior(7,5)*c23;
		P_post(0, 6) = P_prior(0,6)+P_prior(6,6)*c22-P_prior(7,6)*c23;
		P_post(0, 7) = P_prior(0,7)+P_prior(6,7)*c22-P_prior(7,7)*c23;
		P_post(0, 8) = P_prior(0,8)+P_prior(6,8)*c22-P_prior(7,8)*c23;
		P_post(1, 0) = P_post(0, 1);
		P_post(1, 1) = P_prior(1,1)+P_prior(6,1)*c20-P_prior(7,1)*c21;
		P_post(1, 2) = P_prior(1,2)+P_prior(6,2)*c20-P_prior(7,2)*c21;
		P_post(1, 3) = P_prior(1,3)+P_prior(6,3)*c20-P_prior(7,3)*c21;
		P_post(1, 4) = P_prior(1,4)+P_prior(6,4)*c20-P_prior(7,4)*c21;
		P_post(1, 5) = P_prior(1,5)+P_prior(6,5)*c20-P_prior(7,5)*c21;
		P_post(1, 6) = P_prior(1,6)+P_prior(6,6)*c20-P_prior(7,6)*c21;
		P_post(1, 7) = P_prior(1,7)+P_prior(6,7)*c20-P_prior(7,7)*c21;
		P_post(1, 8) = P_prior(1,8)+P_prior(6,8)*c20-P_prior(7,8)*c21;
		P_post(2, 0) = P_post(0, 2);
		P_post(2, 1) = P_post(1, 2);
		P_post(2, 2) = P_prior(2,2)+P_prior(6,2)*c18-P_prior(7,2)*c19;
		P_post(2, 3) = P_prior(2,3)+P_prior(6,3)*c18-P_prior(7,3)*c19;
		P_post(2, 4) = P_prior(2,4)+P_prior(6,4)*c18-P_prior(7,4)*c19;
		P_post(2, 5) = P_prior(2,5)+P_prior(6,5)*c18-P_prior(7,5)*c19;
		P_post(2, 6) = P_prior(2,6)+P_prior(6,6)*c18-P_prior(7,6)*c19;
		P_post(2, 7) = P_prior(2,7)+P_prior(6,7)*c18-P_prior(7,7)*c19;
		P_post(2, 8) = P_prior(2,8)+P_prior(6,8)*c18-P_prior(7,8)*c19;
		P_post(3, 0) = P_post(0, 3);
		P_post(3, 1) = P_post(1, 3);
		P_post(3, 2) = P_post(2, 3);
		P_post(3, 3) = P_prior(3,3)+P_prior(6,3)*c16-P_prior(7,3)*c17;
		P_post(3, 4) = P_prior(3,4)+P_prior(6,4)*c16-P_prior(7,4)*c17;
		P_post(3, 5) = P_prior(3,5)+P_prior(6,5)*c16-P_prior(7,5)*c17;
		P_post(3, 6) = P_prior(3,6)+P_prior(6,6)*c16-P_prior(7,6)*c17;
		P_post(3, 7) = P_prior(3,7)+P_prior(6,7)*c16-P_prior(7,7)*c17;
		P_post(3, 8) = P_prior(3,8)+P_prior(6,8)*c16-P_prior(7,8)*c17;
		P_post(4, 0) = P_post(0, 4);
		P_post(4, 1) = P_post(1, 4);
		P_post(4, 2) = P_post(2, 4);
		P_post(4, 3) = P_post(3, 4);
		P_post(4, 4) = P_prior(4,4)+P_prior(6,4)*c14-P_prior(7,4)*c15;
		P_post(4, 5) = P_prior(4,5)+P_prior(6,5)*c14-P_prior(7,5)*c15;
		P_post(4, 6) = P_prior(4,6)+P_prior(6,6)*c14-P_prior(7,6)*c15;
		P_post(4, 7) = P_prior(4,7)+P_prior(6,7)*c14-P_prior(7,7)*c15;
		P_post(4, 8) = P_prior(4,8)+P_prior(6,8)*c14-P_prior(7,8)*c15;
		P_post(5, 0) = P_post(0, 5);
		P_post(5, 1) = P_post(1, 5);
		P_post(5, 2) = P_post(2, 5);
		P_post(5, 3) = P_post(3, 5);
		P_post(5, 4) = P_post(4, 5);
		P_post(5, 5) = P_prior(5,5)+P_prior(6,5)*c12-P_prior(7,5)*c13;
		P_post(5, 6) = P_prior(5,6)+P_prior(6,6)*c12-P_prior(7,6)*c13;
		P_post(5, 7) = P_prior(5,7)+P_prior(6,7)*c12-P_prior(7,7)*c13;
		P_post(5, 8) = P_prior(5,8)+P_prior(6,8)*c12-P_prior(7,8)*c13;
		P_post(6, 0) = P_post(0, 6);
		P_post(6, 1) = P_post(1, 6);
		P_post(6, 2) = P_post(2, 6);
		P_post(6, 3) = P_post(3, 6);
		P_post(6, 4) = P_post(4, 6);
		P_post(6, 5) = P_post(5, 6);
		P_post(6, 6) = P_prior(6,6)*c7-P_prior(7,6)*c11;
		P_post(6, 7) = P_prior(6,7)*c7-P_prior(7,7)*c11;
		P_post(6, 8) = P_prior(6,8)*c7-P_prior(7,8)*c11;
		P_post(7, 0) = P_post(0, 7);
		P_post(7, 1) = P_post(1, 7);
		P_post(7, 2) = P_post(2, 7);
		P_post(7, 3) = P_post(3, 7);
		P_post(7, 4) = P_post(4, 7);
		P_post(7, 5) = P_post(5, 7);
		P_post(7, 6) = P_post(6, 7);
		P_post(7, 7) = P_prior(6,7)*c10-P_prior(7,7)*c6;
		P_post(7, 8) = P_prior(6,8)*c10-P_prior(7,8)*c6;
		P_post(8, 0) = P_post(0, 8);
		P_post(8, 1) = P_post(1, 8);
		P_post(8, 2) = P_post(2, 8);
		P_post(8, 3) = P_post(3, 8);
		P_post(8, 4) = P_post(4, 8);
		P_post(8, 5) = P_post(5, 8);
		P_post(8, 6) = P_post(6, 8);
		P_post(8, 7) = P_post(7, 8);
		P_post(8, 8) = P_prior(8,8)+P_prior(6,8)*c8-P_prior(7,8)*c9;
	}

	/* _P_post becoms the _P_prior of the magnatometer correction */
	memcpy(mat_data(_P_prior), mat_data(_P_post), sizeof(float) * 9);

	/* error state injection */
	float q_error[4];
	q_error[0] = 1.0f;
	q_error[1] = 0.5 * mat_data(x_error_state)[0];
	q_error[2] = 0.5 * mat_data(x_error_state)[1];
	q_error[3] = 0.5 * mat_data(x_error_state)[2];

	//x_nominal (a posteriori) = q_error * x_nominal (a priori)
	float x_last[4];
	quaternion_copy(x_last, mat_data(x_nominal));
	quaternion_mult(x_last, q_error, mat_data(x_nominal));

	//renormailization
	quat_normalize(mat_data(x_nominal));
}

void eskf_ins_magnetometer_correct(float *mag)
{
	float mx = mag[0];
	float my = mag[1];
	float mz = mag[2];

	float q0 = mat_data(x_nominal)[0];
	float q1 = mat_data(x_nominal)[1];
	float q2 = mat_data(x_nominal)[2];
	float q3 = mat_data(x_nominal)[3];

	float gamma = sqrt(mag[0]*mag[0] + mag[1]*mag[1]);

	/* codeblock for preventing nameing conflict */
	{
		float c0 = P_post(5,8)+P_post(6,8)*R_am_ab_dt(2,0)+P_post(7,8)*R_am_ab_dt(2,1)+P_post(8,8)*R_am_ab_dt(2,2);
		float c1 = P_post(5,7)+P_post(6,7)*R_am_ab_dt(2,0)+P_post(7,7)*R_am_ab_dt(2,1)+P_post(8,7)*R_am_ab_dt(2,2);
		float c2 = P_post(5,6)+P_post(6,6)*R_am_ab_dt(2,0)+P_post(7,6)*R_am_ab_dt(2,1)+P_post(8,6)*R_am_ab_dt(2,2);
		float c3 = P_post(4,8)+P_post(6,8)*R_am_ab_dt(1,0)+P_post(7,8)*R_am_ab_dt(1,1)+P_post(8,8)*R_am_ab_dt(1,2);
		float c4 = P_post(4,7)+P_post(6,7)*R_am_ab_dt(1,0)+P_post(7,7)*R_am_ab_dt(1,1)+P_post(8,7)*R_am_ab_dt(1,2);
		float c5 = P_post(4,6)+P_post(6,6)*R_am_ab_dt(1,0)+P_post(7,6)*R_am_ab_dt(1,1)+P_post(8,6)*R_am_ab_dt(1,2);
		float c6 = P_post(3,8)+P_post(6,8)*R_am_ab_dt(0,0)+P_post(7,8)*R_am_ab_dt(0,1)+P_post(8,8)*R_am_ab_dt(0,2);
		float c7 = P_post(3,7)+P_post(6,7)*R_am_ab_dt(0,0)+P_post(7,7)*R_am_ab_dt(0,1)+P_post(8,7)*R_am_ab_dt(0,2);
		float c8 = P_post(3,6)+P_post(6,6)*R_am_ab_dt(0,0)+P_post(7,6)*R_am_ab_dt(0,1)+P_post(8,6)*R_am_ab_dt(0,2);
		float c9 = P_post(6,8)*Rt_wm_wb_dt(2,0)+P_post(7,8)*Rt_wm_wb_dt(2,1)+P_post(8,8)*Rt_wm_wb_dt(2,2);
		float c10 = P_post(6,7)*Rt_wm_wb_dt(2,0)+P_post(7,7)*Rt_wm_wb_dt(2,1)+P_post(8,7)*Rt_wm_wb_dt(2,2);
		float c11 = P_post(6,6)*Rt_wm_wb_dt(2,0)+P_post(7,6)*Rt_wm_wb_dt(2,1)+P_post(8,6)*Rt_wm_wb_dt(2,2);
		float c12 = P_post(6,8)*Rt_wm_wb_dt(1,0)+P_post(7,8)*Rt_wm_wb_dt(1,1)+P_post(8,8)*Rt_wm_wb_dt(1,2);
		float c13 = P_post(6,7)*Rt_wm_wb_dt(1,0)+P_post(7,7)*Rt_wm_wb_dt(1,1)+P_post(8,7)*Rt_wm_wb_dt(1,2);
		float c14 = P_post(6,6)*Rt_wm_wb_dt(1,0)+P_post(7,6)*Rt_wm_wb_dt(1,1)+P_post(8,6)*Rt_wm_wb_dt(1,2);
		float c15 = P_post(6,8)*Rt_wm_wb_dt(0,0)+P_post(7,8)*Rt_wm_wb_dt(0,1)+P_post(8,8)*Rt_wm_wb_dt(0,2);
		float c16 = P_post(6,7)*Rt_wm_wb_dt(0,0)+P_post(7,7)*Rt_wm_wb_dt(0,1)+P_post(8,7)*Rt_wm_wb_dt(0,2);
		float c17 = P_post(6,6)*Rt_wm_wb_dt(0,0)+P_post(7,6)*Rt_wm_wb_dt(0,1)+P_post(8,6)*Rt_wm_wb_dt(0,2);
		float c18 = P_post(2,8)+P_post(5,8)*dt;
		float c19 = P_post(2,7)+P_post(5,7)*dt;
		float c20 = P_post(2,6)+P_post(5,6)*dt;
		float c21 = P_post(1,8)+P_post(4,8)*dt;
		float c22 = P_post(1,7)+P_post(4,7)*dt;
		float c23 = P_post(1,6)+P_post(4,6)*dt;
		float c24 = P_post(0,8)+P_post(3,8)*dt;
		float c25 = P_post(0,7)+P_post(3,7)*dt;
		float c26 = P_post(0,6)+P_post(3,6)*dt;
		float c27 = P_post(5,4)*dt;
		float c28 = P_post(5,3)*dt;
		float c29 = P_post(4,4)*dt;
		float c30 = P_post(4,3)*dt;
		float c31 = P_post(3,4)*dt;
		float c32 = P_post(3,3)*dt;
		float c33 = P_post(8,4)*Rt_wm_wb_dt(2,2);
		float c34 = P_post(8,3)*Rt_wm_wb_dt(2,2);
		float c35 = P_post(8,4)*Rt_wm_wb_dt(1,2);
		float c36 = P_post(8,3)*Rt_wm_wb_dt(1,2);
		float c37 = P_post(7,4)*Rt_wm_wb_dt(2,1);
		float c38 = P_post(7,3)*Rt_wm_wb_dt(2,1);
		float c39 = P_post(8,4)*Rt_wm_wb_dt(0,2);
		float c40 = P_post(8,3)*Rt_wm_wb_dt(0,2);
		float c41 = P_post(7,4)*Rt_wm_wb_dt(1,1);
		float c42 = P_post(7,3)*Rt_wm_wb_dt(1,1);
		float c43 = P_post(6,4)*Rt_wm_wb_dt(2,0);
		float c44 = P_post(6,3)*Rt_wm_wb_dt(2,0);
		float c45 = P_post(7,4)*Rt_wm_wb_dt(0,1);
		float c46 = P_post(7,3)*Rt_wm_wb_dt(0,1);
		float c47 = P_post(6,4)*Rt_wm_wb_dt(1,0);
		float c48 = P_post(6,3)*Rt_wm_wb_dt(1,0);
		float c49 = P_post(6,4)*Rt_wm_wb_dt(0,0);
		float c50 = P_post(6,3)*Rt_wm_wb_dt(0,0);
		float c51 = P_post(8,4)*R_am_ab_dt(2,2);
		float c52 = P_post(8,3)*R_am_ab_dt(2,2);
		float c53 = P_post(8,4)*R_am_ab_dt(1,2);
		float c54 = P_post(8,3)*R_am_ab_dt(1,2);
		float c55 = P_post(7,4)*R_am_ab_dt(2,1);
		float c56 = P_post(7,3)*R_am_ab_dt(2,1);
		float c57 = P_post(8,4)*R_am_ab_dt(0,2);
		float c58 = P_post(8,3)*R_am_ab_dt(0,2);
		float c59 = P_post(7,4)*R_am_ab_dt(1,1);
		float c60 = P_post(7,3)*R_am_ab_dt(1,1);
		float c61 = P_post(6,4)*R_am_ab_dt(2,0);
		float c62 = P_post(6,3)*R_am_ab_dt(2,0);
		float c63 = P_post(7,4)*R_am_ab_dt(0,1);
		float c64 = P_post(7,3)*R_am_ab_dt(0,1);
		float c65 = P_post(6,4)*R_am_ab_dt(1,0);
		float c66 = P_post(6,3)*R_am_ab_dt(1,0);
		float c67 = P_post(6,4)*R_am_ab_dt(0,0);
		float c68 = P_post(6,3)*R_am_ab_dt(0,0);
		float c69 = P_post(5,3)+c52+c56+c62;
		float c70 = P_post(4,3)+c54+c60+c66;
		float c71 = P_post(3,3)+c58+c64+c68;
		float c72 = P_post(5,4)+c51+c55+c61;
		float c73 = P_post(4,4)+c53+c59+c65;
		float c74 = P_post(3,4)+c57+c63+c67;
		float c75 = c40+c46+c50;
		float c76 = c39+c45+c49;
		float c77 = c36+c42+c48;
		float c78 = c35+c41+c47;
		float c79 = c34+c38+c44;
		float c80 = c33+c37+c43;
		float c81 = P_post(2,4)+c27;
		float c82 = P_post(2,3)+c28;
		float c83 = P_post(1,4)+c29;
		float c84 = P_post(1,3)+c30;
		float c85 = P_post(0,4)+c31;
		float c86 = P_post(0,3)+c32;

		PHt_mag(0, 0) = P_post(0,0)+P_post(3,0)*dt+c86*dt;
		PHt_mag(0, 1) = P_post(0,1)+P_post(3,1)*dt+c85*dt;
		PHt_mag(0, 2) = c86+R_am_ab_dt(0,0)*c26+R_am_ab_dt(0,1)*c25+R_am_ab_dt(0,2)*c24;
		PHt_mag(0, 3) = c85+R_am_ab_dt(1,0)*c26+R_am_ab_dt(1,1)*c25+R_am_ab_dt(1,2)*c24;
		PHt_mag(1, 0) = P_post(1,0)+P_post(4,0)*dt+c84*dt;
		PHt_mag(1, 1) = P_post(1,1)+P_post(4,1)*dt+c83*dt;
		PHt_mag(1, 2) = c84+R_am_ab_dt(0,0)*c23+R_am_ab_dt(0,1)*c22+R_am_ab_dt(0,2)*c21;
		PHt_mag(1, 3) = c83+R_am_ab_dt(1,0)*c23+R_am_ab_dt(1,1)*c22+R_am_ab_dt(1,2)*c21;
		PHt_mag(2, 0) = P_post(2,0)+P_post(5,0)*dt+c82*dt;
		PHt_mag(2, 1) = P_post(2,1)+P_post(5,1)*dt+c81*dt;
		PHt_mag(2, 2) = c82+R_am_ab_dt(0,0)*c20+R_am_ab_dt(0,1)*c19+R_am_ab_dt(0,2)*c18;
		PHt_mag(2, 3) = c81+R_am_ab_dt(1,0)*c20+R_am_ab_dt(1,1)*c19+R_am_ab_dt(1,2)*c18;
		PHt_mag(3, 0) = P_post(3,0)+P_post(6,0)*R_am_ab_dt(0,0)+P_post(7,0)*R_am_ab_dt(0,1)+P_post(8,0)*R_am_ab_dt(0,2)+c71*dt;
		PHt_mag(3, 1) = P_post(3,1)+P_post(6,1)*R_am_ab_dt(0,0)+P_post(7,1)*R_am_ab_dt(0,1)+P_post(8,1)*R_am_ab_dt(0,2)+c74*dt;
		PHt_mag(3, 2) = Q_i(0,0)+c71+R_am_ab_dt(0,0)*c8+R_am_ab_dt(0,1)*c7+R_am_ab_dt(0,2)*c6;
		PHt_mag(3, 3) = c74+R_am_ab_dt(1,0)*c8+R_am_ab_dt(1,1)*c7+R_am_ab_dt(1,2)*c6;
		PHt_mag(4, 0) = P_post(4,0)+P_post(6,0)*R_am_ab_dt(1,0)+P_post(7,0)*R_am_ab_dt(1,1)+P_post(8,0)*R_am_ab_dt(1,2)+c70*dt;
		PHt_mag(4, 1) = P_post(4,1)+P_post(6,1)*R_am_ab_dt(1,0)+P_post(7,1)*R_am_ab_dt(1,1)+P_post(8,1)*R_am_ab_dt(1,2)+c73*dt;
		PHt_mag(4, 2) = c70+R_am_ab_dt(0,0)*c5+R_am_ab_dt(0,1)*c4+R_am_ab_dt(0,2)*c3;
		PHt_mag(4, 3) = Q_i(1,1)+c73+R_am_ab_dt(1,0)*c5+R_am_ab_dt(1,1)*c4+R_am_ab_dt(1,2)*c3;
		PHt_mag(5, 0) = P_post(5,0)+P_post(6,0)*R_am_ab_dt(2,0)+P_post(7,0)*R_am_ab_dt(2,1)+P_post(8,0)*R_am_ab_dt(2,2)+c69*dt;
		PHt_mag(5, 1) = P_post(5,1)+P_post(6,1)*R_am_ab_dt(2,0)+P_post(7,1)*R_am_ab_dt(2,1)+P_post(8,1)*R_am_ab_dt(2,2)+c72*dt;
		PHt_mag(5, 2) = c69+R_am_ab_dt(0,0)*c2+R_am_ab_dt(0,1)*c1+R_am_ab_dt(0,2)*c0;
		PHt_mag(5, 3) = c72+R_am_ab_dt(1,0)*c2+R_am_ab_dt(1,1)*c1+R_am_ab_dt(1,2)*c0;
		PHt_mag(6, 0) = P_post(6,0)*Rt_wm_wb_dt(0,0)+P_post(7,0)*Rt_wm_wb_dt(0,1)+P_post(8,0)*Rt_wm_wb_dt(0,2)+c75*dt;
		PHt_mag(6, 1) = P_post(6,1)*Rt_wm_wb_dt(0,0)+P_post(7,1)*Rt_wm_wb_dt(0,1)+P_post(8,1)*Rt_wm_wb_dt(0,2)+c76*dt;
		PHt_mag(6, 2) = c75+R_am_ab_dt(0,0)*c17+R_am_ab_dt(0,1)*c16+R_am_ab_dt(0,2)*c15;
		PHt_mag(6, 3) = c76+R_am_ab_dt(1,0)*c17+R_am_ab_dt(1,1)*c16+R_am_ab_dt(1,2)*c15;
		PHt_mag(7, 0) = P_post(6,0)*Rt_wm_wb_dt(1,0)+P_post(7,0)*Rt_wm_wb_dt(1,1)+P_post(8,0)*Rt_wm_wb_dt(1,2)+c77*dt;
		PHt_mag(7, 1) = P_post(6,1)*Rt_wm_wb_dt(1,0)+P_post(7,1)*Rt_wm_wb_dt(1,1)+P_post(8,1)*Rt_wm_wb_dt(1,2)+c78*dt;
		PHt_mag(7, 2) = c77+R_am_ab_dt(0,0)*c14+R_am_ab_dt(0,1)*c13+R_am_ab_dt(0,2)*c12;
		PHt_mag(7, 3) = c78+R_am_ab_dt(1,0)*c14+R_am_ab_dt(1,1)*c13+R_am_ab_dt(1,2)*c12;
		PHt_mag(8, 0) = P_post(6,0)*Rt_wm_wb_dt(2,0)+P_post(7,0)*Rt_wm_wb_dt(2,1)+P_post(8,0)*Rt_wm_wb_dt(2,2)+c79*dt;
		PHt_mag(8, 1) = P_post(6,1)*Rt_wm_wb_dt(2,0)+P_post(7,1)*Rt_wm_wb_dt(2,1)+P_post(8,1)*Rt_wm_wb_dt(2,2)+c80*dt;
		PHt_mag(8, 2) = c79+R_am_ab_dt(0,0)*c11+R_am_ab_dt(0,1)*c10+R_am_ab_dt(0,2)*c9;
		PHt_mag(8, 3) = c80+R_am_ab_dt(1,0)*c11+R_am_ab_dt(1,1)*c10+R_am_ab_dt(1,2)*c9;

	}

	/* codeblock for preventing nameing conflict */
	{
	}

	/* codeblock for preventing nameing conflict */
	{
	}

	/* calculate kalman gain */
	//K = P * Ht * inv(H*P*Ht + V)
	MAT_INV(&_HPHt_V_mag, &_HPHt_V_mag_inv);
	MAT_MULT(&_PHt_mag, &_HPHt_V_mag_inv, &_K_mag);

	/* codeblock for preventing nameing conflict */
	{
		/* calculate error state residual */
		float c0_ = -q2*q2;
		float c1_ = q3*q3;
		float c2_ = q1*q1;
		float c3_ = q0*q0;
		float c4_ = q0*q2;
		float c5_ = q1*q3;

		float c0 = mz-mz*(c0_+c1_-c2_+c3_)+gamma*(c4_-c5_)*2.0;
		float c1 = -mx+gamma*(c0_-c1_+c2_+c3_)+mz*(c4_+c5_)*2.0;
		float c2 = my-gamma*(q0*q3+q1*q2)*2.0+mz*(q0*q1-q2*q3)*2.0;

		mat_data(x_error_state)[6] = -K_mag(6,0)*c1+K_mag(6,2)*c0+K_mag(6,1)*c2;
		mat_data(x_error_state)[7] = -K_mag(7,0)*c1+K_mag(7,2)*c0+K_mag(7,1)*c2;
		mat_data(x_error_state)[8] = -K_mag(8,0)*c1+K_mag(8,2)*c0+K_mag(8,1)*c2;
	}

	/* codeblock for preventing nameing conflict */
	{
		/* calculate a posteriori process covariance matrix */
		//P = (I - K*H) * P
		float c0 = gamma*q3*2.0-mz*q1*2.0;
		float c0_ = (c0*q2)*0.5;
		float c7_ = (c0*q3)*0.5;
		float c8_ = (c0*q1)*0.5;
		float c9_ = (c0*q0)*0.5;
		float c1 = gamma*q2*2.0-mz*q0*2.0;
		float c1_ = (c1*q3)*0.5;
		float c5_ = c1*q2*(-1.0*0.5);
		float c14_ = (c1*q1)*0.5;
		float c15_ = (c1*q0)*0.5;
		float c2 = gamma*q1*2.0+mz*q3*2.0;
		float c2_ = (c2*q0)*0.5;
		float c4_ = c2*q2*(-1.0*0.5);
		float c11_ = (c2*q3)*0.5;
		float c13_ = (c2*q1)*0.5;
		float c3 = gamma*q0*2.0+mz*q2*2.0;
		float c3_ = c3*q1*(-1.0*0.5);
		float c6_ = (c3*q2)*0.5;
		float c10_ = (c3*q3)*0.5;
		float c12_ = (c3*q0)*0.5;
		float c4 = -c11_-c15_+c6_+c8_;
		float c5 = -c12_+c13_+c5_+c7_;
		float c6 = c12_-c13_+c5_+c7_;
		float c7 = c0_-c1_+c2_+c3_;
		float c8 = c0_+c1_-c2_+c3_;
		float c9 = -c10_+c14_+c4_+c9_;
		float c10 = c10_-c14_+c4_+c9_;
		float c11 = c0_+c1_+c2_+(c3*q1)*0.5;
		float c12 = c11_+c15_+c6_+c8_;
		float c13 = K_mag(6,1)*c4-K_mag(6,0)*c7-K_mag(6,2)*c9+1.0;
		float c14 = K_mag(8,1)*c5+K_mag(8,0)*c10-K_mag(8,2)*c11+1.0;
		float c15 = -K_mag(8,1)*c4+K_mag(8,0)*c7+K_mag(8,2)*c9;
		float c16 = K_mag(7,1)*c5+K_mag(7,0)*c10-K_mag(7,2)*c11;
		float c17 = -K_mag(7,1)*c4+K_mag(7,0)*c7+K_mag(7,2)*c9;
		float c18 = K_mag(6,1)*c5+K_mag(6,0)*c10-K_mag(6,2)*c11;
		float c19 = K_mag(5,1)*c5+K_mag(5,0)*c10-K_mag(5,2)*c11;
		float c20 = -K_mag(5,1)*c4+K_mag(5,0)*c7+K_mag(5,2)*c9;
		float c21 = K_mag(4,1)*c5+K_mag(4,0)*c10-K_mag(4,2)*c11;
		float c22 = -K_mag(4,1)*c4+K_mag(4,0)*c7+K_mag(4,2)*c9;
		float c23 = K_mag(3,1)*c5+K_mag(3,0)*c10-K_mag(3,2)*c11;
		float c24 = -K_mag(3,1)*c4+K_mag(3,0)*c7+K_mag(3,2)*c9;
		float c25 = K_mag(2,1)*c5+K_mag(2,0)*c10-K_mag(2,2)*c11;
		float c26 = -K_mag(2,1)*c4+K_mag(2,0)*c7+K_mag(2,2)*c9;
		float c27 = K_mag(1,1)*c5+K_mag(1,0)*c10-K_mag(1,2)*c11;
		float c28 = -K_mag(1,1)*c4+K_mag(1,0)*c7+K_mag(1,2)*c9;
		float c29 = K_mag(0,1)*c5+K_mag(0,0)*c10-K_mag(0,2)*c11;
		float c30 = -K_mag(0,1)*c4+K_mag(0,0)*c7+K_mag(0,2)*c9;
		float c31 = K_mag(7,2)*c6+K_mag(7,1)*c8+K_mag(7,0)*c12+1.0;
		float c32 = K_mag(8,2)*c6+K_mag(8,1)*c8+K_mag(8,0)*c12;
		float c33 = K_mag(6,2)*c6+K_mag(6,1)*c8+K_mag(6,0)*c12;
		float c34 = K_mag(5,2)*c6+K_mag(5,1)*c8+K_mag(5,0)*c12;
		float c35 = K_mag(4,2)*c6+K_mag(4,1)*c8+K_mag(4,0)*c12;
		float c36 = K_mag(3,2)*c6+K_mag(3,1)*c8+K_mag(3,0)*c12;
		float c37 = K_mag(2,2)*c6+K_mag(2,1)*c8+K_mag(2,0)*c12;
		float c38 = K_mag(1,2)*c6+K_mag(1,1)*c8+K_mag(1,0)*c12;
		float c39 = K_mag(0,2)*c6+K_mag(0,1)*c8+K_mag(0,0)*c12;

		P_post(0, 0) = P_prior(0,0)-P_prior(6,0)*c30+P_prior(7,0)*c39+P_prior(8,0)*c29;
		P_post(0, 1) = P_prior(0,1)-P_prior(6,1)*c30+P_prior(7,1)*c39+P_prior(8,1)*c29;
		P_post(0, 2) = P_prior(0,2)-P_prior(6,2)*c30+P_prior(7,2)*c39+P_prior(8,2)*c29;
		P_post(0, 3) = P_prior(0,3)-P_prior(6,3)*c30+P_prior(7,3)*c39+P_prior(8,3)*c29;
		P_post(0, 4) = P_prior(0,4)-P_prior(6,4)*c30+P_prior(7,4)*c39+P_prior(8,4)*c29;
		P_post(0, 5) = P_prior(0,5)-P_prior(6,5)*c30+P_prior(7,5)*c39+P_prior(8,5)*c29;
		P_post(0, 6) = P_prior(0,6)-P_prior(6,6)*c30+P_prior(7,6)*c39+P_prior(8,6)*c29;
		P_post(0, 7) = P_prior(0,7)-P_prior(6,7)*c30+P_prior(7,7)*c39+P_prior(8,7)*c29;
		P_post(0, 8) = P_prior(0,8)-P_prior(6,8)*c30+P_prior(7,8)*c39+P_prior(8,8)*c29;
		P_post(1, 0) = P_post(0, 1);
		P_post(1, 1) = P_prior(1,1)-P_prior(6,1)*c28+P_prior(8,1)*c27+P_prior(7,1)*c38;
		P_post(1, 2) = P_prior(1,2)-P_prior(6,2)*c28+P_prior(8,2)*c27+P_prior(7,2)*c38;
		P_post(1, 3) = P_prior(1,3)-P_prior(6,3)*c28+P_prior(8,3)*c27+P_prior(7,3)*c38;
		P_post(1, 4) = P_prior(1,4)-P_prior(6,4)*c28+P_prior(8,4)*c27+P_prior(7,4)*c38;
		P_post(1, 5) = P_prior(1,5)-P_prior(6,5)*c28+P_prior(8,5)*c27+P_prior(7,5)*c38;
		P_post(1, 6) = P_prior(1,6)-P_prior(6,6)*c28+P_prior(8,6)*c27+P_prior(7,6)*c38;
		P_post(1, 7) = P_prior(1,7)-P_prior(6,7)*c28+P_prior(8,7)*c27+P_prior(7,7)*c38;
		P_post(1, 8) = P_prior(1,8)-P_prior(6,8)*c28+P_prior(8,8)*c27+P_prior(7,8)*c38;
		P_post(2, 0) = P_post(0, 2);
		P_post(2, 1) = P_post(1, 2);
		P_post(2, 2) = P_prior(2,2)-P_prior(6,2)*c26+P_prior(8,2)*c25+P_prior(7,2)*c37;
		P_post(2, 3) = P_prior(2,3)-P_prior(6,3)*c26+P_prior(8,3)*c25+P_prior(7,3)*c37;
		P_post(2, 4) = P_prior(2,4)-P_prior(6,4)*c26+P_prior(8,4)*c25+P_prior(7,4)*c37;
		P_post(2, 5) = P_prior(2,5)-P_prior(6,5)*c26+P_prior(8,5)*c25+P_prior(7,5)*c37;
		P_post(2, 6) = P_prior(2,6)-P_prior(6,6)*c26+P_prior(8,6)*c25+P_prior(7,6)*c37;
		P_post(2, 7) = P_prior(2,7)-P_prior(6,7)*c26+P_prior(8,7)*c25+P_prior(7,7)*c37;
		P_post(2, 8) = P_prior(2,8)-P_prior(6,8)*c26+P_prior(8,8)*c25+P_prior(7,8)*c37;
		P_post(3, 0) = P_post(0, 3);
		P_post(3, 1) = P_post(1, 3);
		P_post(3, 2) = P_post(2, 3);
		P_post(3, 3) = P_prior(3,3)-P_prior(6,3)*c24+P_prior(8,3)*c23+P_prior(7,3)*c36;
		P_post(3, 4) = P_prior(3,4)-P_prior(6,4)*c24+P_prior(8,4)*c23+P_prior(7,4)*c36;
		P_post(3, 5) = P_prior(3,5)-P_prior(6,5)*c24+P_prior(8,5)*c23+P_prior(7,5)*c36;
		P_post(3, 6) = P_prior(3,6)-P_prior(6,6)*c24+P_prior(8,6)*c23+P_prior(7,6)*c36;
		P_post(3, 7) = P_prior(3,7)-P_prior(6,7)*c24+P_prior(8,7)*c23+P_prior(7,7)*c36;
		P_post(3, 8) = P_prior(3,8)-P_prior(6,8)*c24+P_prior(8,8)*c23+P_prior(7,8)*c36;
		P_post(4, 0) = P_post(0, 4);
		P_post(4, 1) = P_post(1, 4);
		P_post(4, 2) = P_post(2, 4);
		P_post(4, 3) = P_post(3, 4);
		P_post(4, 4) = P_prior(4,4)-P_prior(6,4)*c22+P_prior(8,4)*c21+P_prior(7,4)*c35;
		P_post(4, 5) = P_prior(4,5)-P_prior(6,5)*c22+P_prior(8,5)*c21+P_prior(7,5)*c35;
		P_post(4, 6) = P_prior(4,6)-P_prior(6,6)*c22+P_prior(8,6)*c21+P_prior(7,6)*c35;
		P_post(4, 7) = P_prior(4,7)-P_prior(6,7)*c22+P_prior(8,7)*c21+P_prior(7,7)*c35;
		P_post(4, 8) = P_prior(4,8)-P_prior(6,8)*c22+P_prior(8,8)*c21+P_prior(7,8)*c35;
		P_post(5, 0) = P_post(0, 5);
		P_post(5, 1) = P_post(1, 5);
		P_post(5, 2) = P_post(2, 5);
		P_post(5, 3) = P_post(3, 5);
		P_post(5, 4) = P_post(4, 5);
		P_post(5, 5) = P_prior(5,5)-P_prior(6,5)*c20+P_prior(8,5)*c19+P_prior(7,5)*c34;
		P_post(5, 6) = P_prior(5,6)-P_prior(6,6)*c20+P_prior(8,6)*c19+P_prior(7,6)*c34;
		P_post(5, 7) = P_prior(5,7)-P_prior(6,7)*c20+P_prior(8,7)*c19+P_prior(7,7)*c34;
		P_post(5, 8) = P_prior(5,8)-P_prior(6,8)*c20+P_prior(8,8)*c19+P_prior(7,8)*c34;
		P_post(6, 0) = P_post(0, 6);
		P_post(6, 1) = P_post(1, 6);
		P_post(6, 2) = P_post(2, 6);
		P_post(6, 3) = P_post(3, 6);
		P_post(6, 4) = P_post(4, 6);
		P_post(6, 5) = P_post(5, 6);
		P_post(6, 6) = P_prior(6,6)*c13+P_prior(8,6)*c18+P_prior(7,6)*c33;
		P_post(6, 7) = P_prior(6,7)*c13+P_prior(8,7)*c18+P_prior(7,7)*c33;
		P_post(6, 8) = P_prior(6,8)*c13+P_prior(8,8)*c18+P_prior(7,8)*c33;
		P_post(7, 0) = P_post(0, 7);
		P_post(7, 1) = P_post(1, 7);
		P_post(7, 2) = P_post(2, 7);
		P_post(7, 3) = P_post(3, 7);
		P_post(7, 4) = P_post(4, 7);
		P_post(7, 5) = P_post(5, 7);
		P_post(7, 6) = P_post(6, 7);
		P_post(7, 7) = -P_prior(6,7)*c17+P_prior(8,7)*c16+P_prior(7,7)*c31;
		P_post(7, 8) = -P_prior(6,8)*c17+P_prior(8,8)*c16+P_prior(7,8)*c31;
		P_post(8, 0) = P_post(0, 8);
		P_post(8, 1) = P_post(1, 8);
		P_post(8, 2) = P_post(2, 8);
		P_post(8, 3) = P_post(3, 8);
		P_post(8, 4) = P_post(4, 8);
		P_post(8, 5) = P_post(5, 8);
		P_post(8, 6) = P_post(6, 8);
		P_post(8, 7) = P_post(7, 8);
		P_post(8, 8) = -P_prior(6,8)*c15+P_prior(8,8)*c14+P_prior(7,8)*c32;
	}

	/* error state injection */
	float q_error[4];
	q_error[0] = 1.0f;
	q_error[1] = 0.0f; //0.5 * mat_data(x_error_state)[0];
	q_error[2] = 0.0f; //0.5 * mat_data(x_error_state)[1];
	q_error[3] = 0.5 * mat_data(x_error_state)[2];

	//x_nominal (a posteriori) = q_error * x_nominal (a priori)
	float x_last[4];
	quaternion_copy(x_last, mat_data(x_nominal));
	quaternion_mult(x_last, q_error, mat_data(x_nominal));

	//renormailization
	quat_normalize(mat_data(x_nominal));
}

void eskf_ins_gps_correct(float px_enu, float py_enu,
                          float vx_enu, float vy_enu)
{
	float px_gps = px_enu;
	float py_gps = py_enu;
	float vx_gps = vx_enu;
	float vy_gps = vy_enu;
	float px = mat_data(x_nominal)[0];
	float py = mat_data(x_nominal)[1];
	float vx = mat_data(x_nominal)[3];
	float vy = mat_data(x_nominal)[4];

	/* codeblock for preventing nameing conflict */
	{
	}

	/* codeblock for preventing nameing conflict */
	{
	}

	/* calculate kalman gain */
	//K = P * Ht * inv(H*P*Ht + V)
	MAT_INV(&_HPHt_V_gps, &_HPHt_V_gps_inv);
	MAT_MULT(&_PHt_gps, &_HPHt_V_gps_inv, &_K_gps);

	/* codeblock for preventing nameing conflict */
	{
		/* calculate error state residual */
		float c0 = vy-vy_gps;
		float c1 = vx-vx_gps;
		float c2 = py-py_gps;
		float c3 = px-px_gps;

		mat_data(x_error_state)[0] = -K_gps(0,0)*c3-K_gps(0,1)*c2-K_gps(0,2)*c1-K_gps(0,3)*c0;
		mat_data(x_error_state)[1] = -K_gps(1,0)*c3-K_gps(1,1)*c2-K_gps(1,2)*c1-K_gps(1,3)*c0;
		mat_data(x_error_state)[3] = -K_gps(3,0)*c3-K_gps(3,1)*c2-K_gps(3,2)*c1-K_gps(3,3)*c0;
		mat_data(x_error_state)[4] = -K_gps(4,0)*c3-K_gps(4,1)*c2-K_gps(4,2)*c1-K_gps(4,3)*c0;
	}

	/* codeblock for preventing nameing conflict */
	{
		/* calculate a posteriori process covariance matrix */
		float c0 = K_gps(4,3)-1.0;
		float c1 = K_gps(3,2)-1.0;
		float c2 = K_gps(1,1)-1.0;
		float c3 = K_gps(0,0)-1.0;

		P_post(0, 0) = -P_prior(0,0)*c3-K_gps(0,1)*P_prior(1,0)-K_gps(0,2)*P_prior(3,0)-K_gps(0,3)*P_prior(4,0);
		P_post(0, 1) = -P_prior(0,1)*c3-K_gps(0,1)*P_prior(1,1)-K_gps(0,2)*P_prior(3,1)-K_gps(0,3)*P_prior(4,1);
		P_post(0, 2) = -P_prior(0,2)*c3-K_gps(0,1)*P_prior(1,2)-K_gps(0,2)*P_prior(3,2)-K_gps(0,3)*P_prior(4,2);
		P_post(0, 3) = -P_prior(0,3)*c3-K_gps(0,1)*P_prior(1,3)-K_gps(0,2)*P_prior(3,3)-K_gps(0,3)*P_prior(4,3);
		P_post(0, 4) = -P_prior(0,4)*c3-K_gps(0,1)*P_prior(1,4)-K_gps(0,2)*P_prior(3,4)-K_gps(0,3)*P_prior(4,4);
		P_post(0, 5) = -P_prior(0,5)*c3-K_gps(0,1)*P_prior(1,5)-K_gps(0,2)*P_prior(3,5)-K_gps(0,3)*P_prior(4,5);
		P_post(0, 6) = -P_prior(0,6)*c3-K_gps(0,1)*P_prior(1,6)-K_gps(0,2)*P_prior(3,6)-K_gps(0,3)*P_prior(4,6);
		P_post(0, 7) = -P_prior(0,7)*c3-K_gps(0,1)*P_prior(1,7)-K_gps(0,2)*P_prior(3,7)-K_gps(0,3)*P_prior(4,7);
		P_post(0, 8) = -P_prior(0,8)*c3-K_gps(0,1)*P_prior(1,8)-K_gps(0,2)*P_prior(3,8)-K_gps(0,3)*P_prior(4,8);
		P_post(1, 0) = P_post(0, 1);
		P_post(1, 1) = -P_prior(1,1)*c2-K_gps(1,0)*P_prior(0,1)-K_gps(1,2)*P_prior(3,1)-K_gps(1,3)*P_prior(4,1);
		P_post(1, 2) = -P_prior(1,2)*c2-K_gps(1,0)*P_prior(0,2)-K_gps(1,2)*P_prior(3,2)-K_gps(1,3)*P_prior(4,2);
		P_post(1, 3) = -P_prior(1,3)*c2-K_gps(1,0)*P_prior(0,3)-K_gps(1,2)*P_prior(3,3)-K_gps(1,3)*P_prior(4,3);
		P_post(1, 4) = -P_prior(1,4)*c2-K_gps(1,0)*P_prior(0,4)-K_gps(1,2)*P_prior(3,4)-K_gps(1,3)*P_prior(4,4);
		P_post(1, 5) = -P_prior(1,5)*c2-K_gps(1,0)*P_prior(0,5)-K_gps(1,2)*P_prior(3,5)-K_gps(1,3)*P_prior(4,5);
		P_post(1, 6) = -P_prior(1,6)*c2-K_gps(1,0)*P_prior(0,6)-K_gps(1,2)*P_prior(3,6)-K_gps(1,3)*P_prior(4,6);
		P_post(1, 7) = -P_prior(1,7)*c2-K_gps(1,0)*P_prior(0,7)-K_gps(1,2)*P_prior(3,7)-K_gps(1,3)*P_prior(4,7);
		P_post(1, 8) = -P_prior(1,8)*c2-K_gps(1,0)*P_prior(0,8)-K_gps(1,2)*P_prior(3,8)-K_gps(1,3)*P_prior(4,8);
		P_post(2, 0) = P_post(0, 2);
		P_post(2, 1) = P_post(1, 2);
		P_post(2, 2) = P_prior(2,2)-K_gps(2,0)*P_prior(0,2)-K_gps(2,1)*P_prior(1,2)-K_gps(2,2)*P_prior(3,2)-K_gps(2,3)*P_prior(4,2);
		P_post(2, 3) = P_prior(2,3)-K_gps(2,0)*P_prior(0,3)-K_gps(2,1)*P_prior(1,3)-K_gps(2,2)*P_prior(3,3)-K_gps(2,3)*P_prior(4,3);
		P_post(2, 4) = P_prior(2,4)-K_gps(2,0)*P_prior(0,4)-K_gps(2,1)*P_prior(1,4)-K_gps(2,2)*P_prior(3,4)-K_gps(2,3)*P_prior(4,4);
		P_post(2, 5) = P_prior(2,5)-K_gps(2,0)*P_prior(0,5)-K_gps(2,1)*P_prior(1,5)-K_gps(2,2)*P_prior(3,5)-K_gps(2,3)*P_prior(4,5);
		P_post(2, 6) = P_prior(2,6)-K_gps(2,0)*P_prior(0,6)-K_gps(2,1)*P_prior(1,6)-K_gps(2,2)*P_prior(3,6)-K_gps(2,3)*P_prior(4,6);
		P_post(2, 7) = P_prior(2,7)-K_gps(2,0)*P_prior(0,7)-K_gps(2,1)*P_prior(1,7)-K_gps(2,2)*P_prior(3,7)-K_gps(2,3)*P_prior(4,7);
		P_post(2, 8) = P_prior(2,8)-K_gps(2,0)*P_prior(0,8)-K_gps(2,1)*P_prior(1,8)-K_gps(2,2)*P_prior(3,8)-K_gps(2,3)*P_prior(4,8);
		P_post(3, 0) = P_post(0, 3);
		P_post(3, 1) = P_post(1, 3);
		P_post(3, 2) = P_post(2, 3);
		P_post(3, 3) = -P_prior(3,3)*c1-K_gps(3,0)*P_prior(0,3)-K_gps(3,1)*P_prior(1,3)-K_gps(3,3)*P_prior(4,3);
		P_post(3, 4) = -P_prior(3,4)*c1-K_gps(3,0)*P_prior(0,4)-K_gps(3,1)*P_prior(1,4)-K_gps(3,3)*P_prior(4,4);
		P_post(3, 5) = -P_prior(3,5)*c1-K_gps(3,0)*P_prior(0,5)-K_gps(3,1)*P_prior(1,5)-K_gps(3,3)*P_prior(4,5);
		P_post(3, 6) = -P_prior(3,6)*c1-K_gps(3,0)*P_prior(0,6)-K_gps(3,1)*P_prior(1,6)-K_gps(3,3)*P_prior(4,6);
		P_post(3, 7) = -P_prior(3,7)*c1-K_gps(3,0)*P_prior(0,7)-K_gps(3,1)*P_prior(1,7)-K_gps(3,3)*P_prior(4,7);
		P_post(3, 8) = -P_prior(3,8)*c1-K_gps(3,0)*P_prior(0,8)-K_gps(3,1)*P_prior(1,8)-K_gps(3,3)*P_prior(4,8);
		P_post(4, 0) = P_post(0, 4);
		P_post(4, 1) = P_post(1, 4);
		P_post(4, 2) = P_post(2, 4);
		P_post(4, 3) = P_post(3, 4);
		P_post(4, 4) = -P_prior(4,4)*c0-K_gps(4,0)*P_prior(0,4)-K_gps(4,1)*P_prior(1,4)-K_gps(4,2)*P_prior(3,4);
		P_post(4, 5) = -P_prior(4,5)*c0-K_gps(4,0)*P_prior(0,5)-K_gps(4,1)*P_prior(1,5)-K_gps(4,2)*P_prior(3,5);
		P_post(4, 6) = -P_prior(4,6)*c0-K_gps(4,0)*P_prior(0,6)-K_gps(4,1)*P_prior(1,6)-K_gps(4,2)*P_prior(3,6);
		P_post(4, 7) = -P_prior(4,7)*c0-K_gps(4,0)*P_prior(0,7)-K_gps(4,1)*P_prior(1,7)-K_gps(4,2)*P_prior(3,7);
		P_post(4, 8) = -P_prior(4,8)*c0-K_gps(4,0)*P_prior(0,8)-K_gps(4,1)*P_prior(1,8)-K_gps(4,2)*P_prior(3,8);
		P_post(5, 0) = P_post(0, 5);
		P_post(5, 1) = P_post(1, 5);
		P_post(5, 2) = P_post(2, 5);
		P_post(5, 3) = P_post(3, 5);
		P_post(5, 4) = P_post(4, 5);
		P_post(5, 5) = P_prior(5,5)-K_gps(5,0)*P_prior(0,5)-K_gps(5,1)*P_prior(1,5)-K_gps(5,2)*P_prior(3,5)-K_gps(5,3)*P_prior(4,5);
		P_post(5, 6) = P_prior(5,6)-K_gps(5,0)*P_prior(0,6)-K_gps(5,1)*P_prior(1,6)-K_gps(5,2)*P_prior(3,6)-K_gps(5,3)*P_prior(4,6);
		P_post(5, 7) = P_prior(5,7)-K_gps(5,0)*P_prior(0,7)-K_gps(5,1)*P_prior(1,7)-K_gps(5,2)*P_prior(3,7)-K_gps(5,3)*P_prior(4,7);
		P_post(5, 8) = P_prior(5,8)-K_gps(5,0)*P_prior(0,8)-K_gps(5,1)*P_prior(1,8)-K_gps(5,2)*P_prior(3,8)-K_gps(5,3)*P_prior(4,8);
		P_post(6, 0) = P_post(0, 6);
		P_post(6, 1) = P_post(1, 6);
		P_post(6, 2) = P_post(2, 6);
		P_post(6, 3) = P_post(3, 6);
		P_post(6, 4) = P_post(4, 6);
		P_post(6, 5) = P_post(5, 6);
		P_post(6, 6) = P_prior(6,6)-K_gps(6,0)*P_prior(0,6)-K_gps(6,1)*P_prior(1,6)-K_gps(6,2)*P_prior(3,6)-K_gps(6,3)*P_prior(4,6);
		P_post(6, 7) = P_prior(6,7)-K_gps(6,0)*P_prior(0,7)-K_gps(6,1)*P_prior(1,7)-K_gps(6,2)*P_prior(3,7)-K_gps(6,3)*P_prior(4,7);
		P_post(6, 8) = P_prior(6,8)-K_gps(6,0)*P_prior(0,8)-K_gps(6,1)*P_prior(1,8)-K_gps(6,2)*P_prior(3,8)-K_gps(6,3)*P_prior(4,8);
		P_post(7, 0) = P_post(0, 7);
		P_post(7, 1) = P_post(1, 7);
		P_post(7, 2) = P_post(2, 7);
		P_post(7, 3) = P_post(3, 7);
		P_post(7, 4) = P_post(4, 7);
		P_post(7, 5) = P_post(5, 7);
		P_post(7, 6) = P_post(6, 7);
		P_post(7, 7) = P_prior(7,7)-K_gps(7,0)*P_prior(0,7)-K_gps(7,1)*P_prior(1,7)-K_gps(7,2)*P_prior(3,7)-K_gps(7,3)*P_prior(4,7);
		P_post(7, 8) = P_prior(7,8)-K_gps(7,0)*P_prior(0,8)-K_gps(7,1)*P_prior(1,8)-K_gps(7,2)*P_prior(3,8)-K_gps(7,3)*P_prior(4,8);
		P_post(8, 0) = P_post(0, 8);
		P_post(8, 1) = P_post(1, 8);
		P_post(8, 2) = P_post(2, 8);
		P_post(8, 3) = P_post(3, 8);
		P_post(8, 4) = P_post(4, 8);
		P_post(8, 5) = P_post(5, 8);
		P_post(8, 6) = P_post(6, 8);
		P_post(8, 7) = P_post(7, 8);
		P_post(8, 8) = P_prior(8,8)-K_gps(8,0)*P_prior(0,8)-K_gps(8,1)*P_prior(1,8)-K_gps(8,2)*P_prior(3,8)-K_gps(8,3)*P_prior(4,8);
	}

	/* error state injection */
	mat_data(x_nominal)[0] += mat_data(x_error_state)[0];
	mat_data(x_nominal)[1] += mat_data(x_error_state)[1];
	mat_data(x_nominal)[3] += mat_data(x_error_state)[3];
	mat_data(x_nominal)[4] += mat_data(x_error_state)[4];
}

void eskf_ins_barometer_correct(float barometer_z, float barometer_vz)
{
	float pz_baro = barometer_z;
	float vz_baro = barometer_vz;
	float pz = mat_data(x_nominal)[2];
	float vz = mat_data(x_nominal)[5];

	/* codeblock for preventing nameing conflict */
	{
	}

	/* codeblock for preventing nameing conflict */
	{
	}

	/* calculate kalman gain */
	//K = P * Ht * inv(H*P*Ht + V)
	MAT_INV(&_HPHt_V_baro, &_HPHt_V_baro_inv);
	MAT_MULT(&_PHt_baro, &_HPHt_V_baro_inv, &_K_baro);

	/* codeblock for preventing nameing conflict */
	{
		/* calculate error state residual */
		float c0 = vz-vz_baro;
		float c1 = pz-pz_baro;

		mat_data(x_error_state)[2] = -K_baro(2,0)*c1-K_baro(2,1)*c0;
		mat_data(x_error_state)[5] = -K_baro(5,0)*c1-K_baro(5,1)*c0;
	}

	/* codeblock for preventing nameing conflict */
	{
		/* calculate a posteriori process covariance matrix */
		float c0 = K_baro(5,1)-1.0;
		float c1 = K_baro(2,0)-1.0;

		P_post(0, 0) = P_prior(0,0)-K_baro(0,0)*P_prior(2,0)-K_baro(0,1)*P_prior(5,0);
		P_post(0, 1) = P_prior(0,1)-K_baro(0,0)*P_prior(2,1)-K_baro(0,1)*P_prior(5,1);
		P_post(0, 2) = P_prior(0,2)-K_baro(0,0)*P_prior(2,2)-K_baro(0,1)*P_prior(5,2);
		P_post(0, 3) = P_prior(0,3)-K_baro(0,0)*P_prior(2,3)-K_baro(0,1)*P_prior(5,3);
		P_post(0, 4) = P_prior(0,4)-K_baro(0,0)*P_prior(2,4)-K_baro(0,1)*P_prior(5,4);
		P_post(0, 5) = P_prior(0,5)-K_baro(0,0)*P_prior(2,5)-K_baro(0,1)*P_prior(5,5);
		P_post(0, 6) = P_prior(0,6)-K_baro(0,0)*P_prior(2,6)-K_baro(0,1)*P_prior(5,6);
		P_post(0, 7) = P_prior(0,7)-K_baro(0,0)*P_prior(2,7)-K_baro(0,1)*P_prior(5,7);
		P_post(0, 8) = P_prior(0,8)-K_baro(0,0)*P_prior(2,8)-K_baro(0,1)*P_prior(5,8);
		P_post(1, 0) = P_post(0, 1);
		P_post(1, 1) = P_prior(1,1)-K_baro(1,0)*P_prior(2,1)-K_baro(1,1)*P_prior(5,1);
		P_post(1, 2) = P_prior(1,2)-K_baro(1,0)*P_prior(2,2)-K_baro(1,1)*P_prior(5,2);
		P_post(1, 3) = P_prior(1,3)-K_baro(1,0)*P_prior(2,3)-K_baro(1,1)*P_prior(5,3);
		P_post(1, 4) = P_prior(1,4)-K_baro(1,0)*P_prior(2,4)-K_baro(1,1)*P_prior(5,4);
		P_post(1, 5) = P_prior(1,5)-K_baro(1,0)*P_prior(2,5)-K_baro(1,1)*P_prior(5,5);
		P_post(1, 6) = P_prior(1,6)-K_baro(1,0)*P_prior(2,6)-K_baro(1,1)*P_prior(5,6);
		P_post(1, 7) = P_prior(1,7)-K_baro(1,0)*P_prior(2,7)-K_baro(1,1)*P_prior(5,7);
		P_post(1, 8) = P_prior(1,8)-K_baro(1,0)*P_prior(2,8)-K_baro(1,1)*P_prior(5,8);
		P_post(2, 0) = P_post(0, 2);
		P_post(2, 1) = P_post(1, 2);
		P_post(2, 2) = -P_prior(2,2)*c1-K_baro(2,1)*P_prior(5,2);
		P_post(2, 3) = -P_prior(2,3)*c1-K_baro(2,1)*P_prior(5,3);
		P_post(2, 4) = -P_prior(2,4)*c1-K_baro(2,1)*P_prior(5,4);
		P_post(2, 5) = -P_prior(2,5)*c1-K_baro(2,1)*P_prior(5,5);
		P_post(2, 6) = -P_prior(2,6)*c1-K_baro(2,1)*P_prior(5,6);
		P_post(2, 7) = -P_prior(2,7)*c1-K_baro(2,1)*P_prior(5,7);
		P_post(2, 8) = -P_prior(2,8)*c1-K_baro(2,1)*P_prior(5,8);
		P_post(3, 0) = P_post(0, 3);
		P_post(3, 1) = P_post(1, 3);
		P_post(3, 2) = P_post(2, 3);
		P_post(3, 3) = P_prior(3,3)-K_baro(3,0)*P_prior(2,3)-K_baro(3,1)*P_prior(5,3);
		P_post(3, 4) = P_prior(3,4)-K_baro(3,0)*P_prior(2,4)-K_baro(3,1)*P_prior(5,4);
		P_post(3, 5) = P_prior(3,5)-K_baro(3,0)*P_prior(2,5)-K_baro(3,1)*P_prior(5,5);
		P_post(3, 6) = P_prior(3,6)-K_baro(3,0)*P_prior(2,6)-K_baro(3,1)*P_prior(5,6);
		P_post(3, 7) = P_prior(3,7)-K_baro(3,0)*P_prior(2,7)-K_baro(3,1)*P_prior(5,7);
		P_post(3, 8) = P_prior(3,8)-K_baro(3,0)*P_prior(2,8)-K_baro(3,1)*P_prior(5,8);
		P_post(4, 0) = P_post(0, 4);
		P_post(4, 1) = P_post(1, 4);
		P_post(4, 2) = P_post(2, 4);
		P_post(4, 3) = P_post(3, 4);
		P_post(4, 4) = P_prior(4,4)-K_baro(4,0)*P_prior(2,4)-K_baro(4,1)*P_prior(5,4);
		P_post(4, 5) = P_prior(4,5)-K_baro(4,0)*P_prior(2,5)-K_baro(4,1)*P_prior(5,5);
		P_post(4, 6) = P_prior(4,6)-K_baro(4,0)*P_prior(2,6)-K_baro(4,1)*P_prior(5,6);
		P_post(4, 7) = P_prior(4,7)-K_baro(4,0)*P_prior(2,7)-K_baro(4,1)*P_prior(5,7);
		P_post(4, 8) = P_prior(4,8)-K_baro(4,0)*P_prior(2,8)-K_baro(4,1)*P_prior(5,8);
		P_post(5, 0) = P_post(0, 5);
		P_post(5, 1) = P_post(1, 5);
		P_post(5, 2) = P_post(2, 5);
		P_post(5, 3) = P_post(3, 5);
		P_post(5, 4) = P_post(4, 5);
		P_post(5, 5) = -P_prior(5,5)*c0-K_baro(5,0)*P_prior(2,5);
		P_post(5, 6) = -P_prior(5,6)*c0-K_baro(5,0)*P_prior(2,6);
		P_post(5, 7) = -P_prior(5,7)*c0-K_baro(5,0)*P_prior(2,7);
		P_post(5, 8) = -P_prior(5,8)*c0-K_baro(5,0)*P_prior(2,8);
		P_post(6, 0) = P_post(0, 6);
		P_post(6, 1) = P_post(1, 6);
		P_post(6, 2) = P_post(2, 6);
		P_post(6, 3) = P_post(3, 6);
		P_post(6, 4) = P_post(4, 6);
		P_post(6, 5) = P_post(5, 6);
		P_post(6, 6) = P_prior(6,6)-K_baro(6,0)*P_prior(2,6)-K_baro(6,1)*P_prior(5,6);
		P_post(6, 7) = P_prior(6,7)-K_baro(6,0)*P_prior(2,7)-K_baro(6,1)*P_prior(5,7);
		P_post(6, 8) = P_prior(6,8)-K_baro(6,0)*P_prior(2,8)-K_baro(6,1)*P_prior(5,8);
		P_post(7, 0) = P_post(0, 7);
		P_post(7, 1) = P_post(1, 7);
		P_post(7, 2) = P_post(2, 7);
		P_post(7, 3) = P_post(3, 7);
		P_post(7, 4) = P_post(4, 7);
		P_post(7, 5) = P_post(5, 7);
		P_post(7, 6) = P_post(6, 7);
		P_post(7, 7) = P_prior(7,7)-K_baro(7,0)*P_prior(2,7)-K_baro(7,1)*P_prior(5,7);
		P_post(7, 8) = P_prior(7,8)-K_baro(7,0)*P_prior(2,8)-K_baro(7,1)*P_prior(5,8);
		P_post(8, 0) = P_post(0, 8);
		P_post(8, 1) = P_post(1, 8);
		P_post(8, 2) = P_post(2, 8);
		P_post(8, 3) = P_post(3, 8);
		P_post(8, 4) = P_post(4, 8);
		P_post(8, 5) = P_post(5, 8);
		P_post(8, 6) = P_post(6, 8);
		P_post(8, 7) = P_post(7, 8);
		P_post(8, 8) = P_prior(8,8)-K_baro(8,0)*P_prior(2,8)-K_baro(8,1)*P_prior(5,8);
	}

	/* error state injection */
	mat_data(x_nominal)[2] += mat_data(x_error_state)[2];
	mat_data(x_nominal)[5] += mat_data(x_error_state)[5];
}

void get_eskf_ins_attitude_quaternion(float *q_out)
{
	/* return the conjugated quaternion since we use opposite convention compared to the paper.
	 * paper: quaternion of earth frame to body-fixed frame
	 * us: quaternion of body-fixed frame to earth frame */
	quaternion_conj(mat_data(x_nominal), q_out);
}