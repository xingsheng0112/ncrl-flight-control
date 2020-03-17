float vel_last[3] = {0.0f, 0.0f, 0.0f};
float vel_predict[3] = {0.0f, 0.0f, 0.0f};

void nav_velocity_predict(float *dcm_b2i, float *imu_acc)
{
	/* convert acceleration from body-frame to inertial-frame */
	float acc_i_frame[3];
	acc_i_frame[0] = dcm_b2i[0*3 + 0]*imu_acc[0] + dcm_b2i[0*3 + 1]*imu_acc[1] + dcm_b2i[0*3 + 2]*imu_acc[2];
	acc_i_frame[1] = dcm_b2i[1*3 + 0]*imu_acc[0] + dcm_b2i[1*3 + 1]*imu_acc[1] + dcm_b2i[1*3 + 2]*imu_acc[2];
	acc_i_frame[2] = dcm_b2i[2*3 + 0]*imu_acc[0] + dcm_b2i[2*3 + 1]*imu_acc[1] + dcm_b2i[2*3 + 2]*imu_acc[2] - 9.8f;
	acc_i_frame[0] *= 100.0f;
	acc_i_frame[1] *= 100.0f;
	acc_i_frame[2] *= -100.0f;

	/* integrate translational acceleration to velocity */
	const float dt = 0.0025;
	vel_predict[0] = vel_last[0] + acc_i_frame[0] * dt;
	vel_predict[1] = vel_last[1] + acc_i_frame[1] * dt;
	vel_predict[2] = vel_last[2] + acc_i_frame[2] * dt;
	vel_last[0] = vel_predict[0];
	vel_last[1] = vel_predict[1];
	vel_last[2] = vel_predict[2];
}

void nav_velocity_correct(float *vel_ref, float *vel_filtered)
{
	/* complementary filter */
	const float a = 0.1f;
	vel_filtered[0] = a * vel_predict[0] + (1.0f - a) * vel_ref[0];
	vel_filtered[1] = a * vel_predict[1] + (1.0f - a) * vel_ref[1];
	vel_filtered[2] = a * vel_predict[2] + (1.0f - a) * vel_ref[2];
	vel_predict[0] = vel_last[0] = vel_predict[0];
	vel_predict[1] = vel_last[1] = vel_predict[1];
	vel_predict[2] = vel_last[2] = vel_predict[2];
}
