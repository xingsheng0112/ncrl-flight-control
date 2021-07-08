#include "optitrack.h"
#include "vins_mono.h"
#include "system_state.h"
#include "proj_config.h"
#include "ms5611.h"
#include "gps.h"
#include "ins.h"
#include "autopilot.h"
#include "ist8310.h"
#include "barometer.h"
#include "rangefinder.h"
#include "vio.h"

sensor_manager_t sensor_manager = {
	.heading_src = SELECT_HEADING_SENSOR,
	.height_src = SELECT_HEIGHT_SENSOR,
	.position_src = SELECT_POSITION_SENSOR
};

bool is_heading_available(void)
{
	switch(sensor_manager.heading_src) {
	case HEADING_FUSION_USE_COMPASS:
		return ist8310_available();
	case HEADING_FUSION_USE_OPTITRACK:
		return optitrack_available();
	case HEADING_FUSION_USE_VINS_MONO:
		return vio_available();
	default:
		return false;
	}
}

bool is_xy_position_available(void)
{
	switch(sensor_manager.position_src) {
	case POSITION_FUSION_USE_OPTITRACK:
		return optitrack_available();
	case POSITION_FUSION_USE_GPS:
		return is_gps_available();
	case POSITION_FUSION_USE_VINS_MONO:
		return vio_available();
	default:
		return false;
	}
}

bool is_height_available(void)
{
	switch(sensor_manager.height_src) {
	case HEIGHT_FUSION_USE_OPTITRACK:
		return optitrack_available();
	case HEIGHT_FUSION_USE_BAROMETER:
		return is_barometer_available();
	case HEIGHT_FUSION_USE_VINS_MONO:
		return vio_available();
	case HEIGHT_FUSION_USE_RANGEFINDER:
		return rangefinder_available();
	default:
		return false;
	}
}

void get_attitude_euler_angles(float *roll, float *pitch, float *yaw)
{
	ins_ahrs_get_attitude_euler_angles(roll, pitch, yaw);
}

void get_attitude_quaternion(float *q)
{
	ins_ahrs_get_attitude_quaternion(q);
}

void get_rotation_matrix_b2i(float **R_b2i)
{
	ins_ahrs_get_rotation_matrix_b2i(R_b2i);
}

void get_rotation_matrix_i2b(float **R_i2b)
{
	ins_ahrs_get_rotation_matrix_i2b(R_i2b);
}

void get_enu_position(float *pos)
{
	/* x-y position */
	switch(sensor_manager.position_src) {
	case POSITION_FUSION_USE_OPTITRACK:
		pos[0] = optitrack_read_pos_x();
		pos[1] = optitrack_read_pos_y();
		break;
	case POSITION_FUSION_USE_VINS_MONO:
		pos[0] = vio_get_position_x();
		pos[1] = vio_get_position_y();
		break;
	case POSITION_FUSION_USE_GPS:
		pos[0] = ins_get_fused_position_x();
		pos[1] = ins_get_fused_position_y();
		break;
	default:
		pos[0] = 0.0f;
		pos[1] = 0.0f;
		break;
	}

	/* z position */
	switch(sensor_manager.height_src) {
	case HEIGHT_FUSION_USE_OPTITRACK:
		pos[2] = optitrack_read_pos_z();
		break;
	case HEIGHT_FUSION_USE_VINS_MONO:
		pos[2] = vio_get_position_z();
		break;
	case HEIGHT_FUSION_USE_BAROMETER:
		pos[2] = ins_get_fused_position_z();
		break;
	case HEIGHT_FUSION_USE_RANGEFINDER:
		pos[2] = ins_get_fused_position_z();
		break;
	default:
		pos[2] = 0.0f;
		break;
	}
}

float get_enu_position_x(void)
{
	switch(sensor_manager.position_src) {
	case POSITION_FUSION_USE_OPTITRACK:
		return optitrack_read_pos_x();
	case POSITION_FUSION_USE_VINS_MONO:
		return vio_get_position_x();
	case POSITION_FUSION_USE_GPS:
		return ins_get_fused_position_x();
	default:
		return 0.0f;
	}
}

float get_enu_position_y(void)
{
	switch(sensor_manager.position_src) {
	case POSITION_FUSION_USE_OPTITRACK:
		return optitrack_read_pos_y();
	case POSITION_FUSION_USE_VINS_MONO:
		return vio_get_position_y();
	case POSITION_FUSION_USE_GPS:
		return ins_get_fused_position_y();
	default:
		return 0.0f;
	}
}

float get_enu_position_z(void)
{
	switch(sensor_manager.height_src) {
	case HEIGHT_FUSION_USE_OPTITRACK:
		return optitrack_read_pos_z();
	case HEIGHT_FUSION_USE_VINS_MONO:
		return vio_get_position_z();
	case HEIGHT_FUSION_USE_BAROMETER:
		return ins_get_fused_position_z();
	case HEIGHT_FUSION_USE_RANGEFINDER:
		return ins_get_fused_position_z();
	default:
		return 0.0f;
	}
}

void get_enu_velocity(float *vel)
{
	/* x-y velocity */
	switch(sensor_manager.position_src) {
	case POSITION_FUSION_USE_OPTITRACK:
		vel[0] = optitrack_read_vel_x();
		vel[1] = optitrack_read_vel_y();
		break;
	case POSITION_FUSION_USE_VINS_MONO:
		vel[0] = vio_get_velocity_x();
		vel[1] = vio_get_velocity_y();
		break;
	case POSITION_FUSION_USE_GPS:
		vel[0] = ins_get_fused_velocity_x();
		vel[1] = ins_get_fused_velocity_y();
		break;
	default:
		vel[0] = 0.0f;
		vel[1] = 0.0f;
		break;
	}

	/* z velocity */
	switch(sensor_manager.height_src) {
	case HEIGHT_FUSION_USE_OPTITRACK:
		vel[2] = optitrack_read_vel_z();
		break;
	case HEIGHT_FUSION_USE_VINS_MONO:
		vel[2] = vio_get_velocity_z();
		break;
	case HEIGHT_FUSION_USE_BAROMETER:
		vel[2] = ins_get_fused_velocity_z();
		break;
	case HEIGHT_FUSION_USE_RANGEFINDER:
		vel[2] = ins_get_fused_velocity_z();
		break;
	default:
		vel[2] = 0.0f;
		break;
	}
}

float get_enu_velocity_x(void)
{
	switch(sensor_manager.position_src) {
	case POSITION_FUSION_USE_OPTITRACK:
		return optitrack_read_vel_x();
	case POSITION_FUSION_USE_VINS_MONO:
		return vio_get_velocity_x();
	case POSITION_FUSION_USE_GPS:
		return ins_get_fused_velocity_x();
	default:
		return 0.0f;
	}
}

float get_enu_velocity_y(void)
{
	switch(sensor_manager.position_src) {
	case POSITION_FUSION_USE_OPTITRACK:
		return optitrack_read_vel_y();
	case POSITION_FUSION_USE_VINS_MONO:
		return vio_get_velocity_y();
	case POSITION_FUSION_USE_GPS:
		return ins_get_fused_velocity_y();
	default:
		return 0.0f;
	}
}

float get_enu_velocity_z(void)
{
	switch(sensor_manager.height_src) {
	case HEIGHT_FUSION_USE_OPTITRACK:
		return optitrack_read_vel_z();
	case HEIGHT_FUSION_USE_VINS_MONO:
		return vio_get_velocity_z();
	case HEIGHT_FUSION_USE_BAROMETER:
		return ins_get_fused_velocity_z();
	case HEIGHT_FUSION_USE_RANGEFINDER:
		return ins_get_fused_velocity_z();
	default:
		return 0.0f;
	}
}

void change_heading_sensor_src(int new_src)
{
	if(sensor_manager.heading_src == new_src) {
		return;
	}

	/* check if the sensor to switch is currently available */
	bool sensor_available = false;

	switch(new_src) {
	case HEADING_FUSION_USE_COMPASS:
		sensor_available = ist8310_available();
		break;
	case HEADING_FUSION_USE_OPTITRACK:
		sensor_available = optitrack_available();
		break;
	case HEADING_FUSION_USE_VINS_MONO:
		sensor_available = vio_available();
		break;
	default:
		sensor_available = false;
		break;
	}

	if(sensor_available == false) {
		return;
	}

	sensor_manager.heading_src = new_src;
}

void change_height_sensor_src(int new_src)
{
	if(sensor_manager.height_src == new_src) {
		return;
	}

	/* check if the sensor to switch is currently available */
	bool sensor_available = false;

	switch(new_src) {
	case HEIGHT_FUSION_USE_OPTITRACK:
		sensor_available = optitrack_available();
		break;
	case HEIGHT_FUSION_USE_BAROMETER:
		sensor_available = is_barometer_available();
		break;
	case HEIGHT_FUSION_USE_VINS_MONO:
		sensor_available = vio_available();
		break;
	case HEIGHT_FUSION_USE_RANGEFINDER:
		sensor_available = rangefinder_available();
		break;
	default:
		sensor_available = false;
		break;
	}

	if(sensor_available == false) {
		return;
	}
}

void change_position_sensor_src(int new_src)
{
	if(sensor_manager.position_src == new_src) {
		return;
	}

	/* check if the sensor to switch is currently available */
	bool sensor_available = false;

	switch(new_src) {
	case POSITION_FUSION_USE_OPTITRACK:
		sensor_available = optitrack_available();
		break;
	case POSITION_FUSION_USE_GPS:
		sensor_available = is_gps_available();
		break;
	case POSITION_FUSION_USE_VINS_MONO:
		sensor_available = vio_available();
		break;
	default:
		sensor_available = false;
		break;
	}

	if(sensor_available == false) {
		return;
	}
}
