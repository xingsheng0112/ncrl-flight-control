#include <stdbool.h>
#include "FreeRTOS.h"
#include "task.h"
#include "semphr.h"
#include "ist8310.h"
#include "lidar_lite.h"
#include "delay.h"
#include "ins_sensor_sync.h"
#include "proj_config.h"
#include "sw_i2c.h"
#include "timer.h"
#include "flash.h"
#include "uart.h"
#include "led.h"
#include "ublox_m8n.h"
#include "optitrack.h"
#include "vins_mono.h"
#include "spi.h"
#include "crc.h"
#include "pwm.h"
#include "exti.h"

void f4_sw_i2c_driver_register_task(const char *task_name, configSTACK_DEPTH_TYPE stack_size,
                                    UBaseType_t priority);

void f4_board_init(void)
{
	/* driver initialization */
	flash_init();
	_crc_init();
	led_init();
	ext_switch_init();
	uart1_init(115200);
	uart3_init(115200); //telem
	uart4_init(100000); //s-bus

#if (SELECT_NAVIGATION_DEVICE1 == NAV_DEV1_USE_GPS)
	uart7_init(38400); //gps
	ublox_m8n_init();
#elif (SELECT_NAVIGATION_DEVICE1 == NAV_DEV1_USE_OPTITRACK)
	uart7_init(115200);
	optitrack_init(UAV_DEFAULT_ID); //setup tracker id for this MAV
#endif

#if (SELECT_NAVIGATION_DEVICE2 == NAV_DEV2_USE_VINS_MONO)
	uart6_init(115200);
	vins_mono_init(UAV_DEFAULT_ID); //TODO: tracker id is not needed
#endif

	timer12_init();    //system timer and flight controller timer
	pwm_timer1_init(); //motor
	pwm_timer4_init(); //motor
	exti10_init();     //imu ext interrupt
	spi1_init();       //imu

	blocked_delay_ms(50);

#if ((ENABLE_MAGNETOMETER != 0) || (ENABLE_RANGEFINDER != 0))
	sw_i2c_init();
	f4_sw_i2c_driver_register_task("sw i2c driver", 512, tskIDLE_PRIORITY + 5);
#endif

#if (ENABLE_BAROMETER != 0)
	/* barometer (ms5611) */
	spi3_init();
	ms5611_init();
#endif

	/* sensor scheduler timer */
	timer3_init();
}

void f4_sw_i2c_driver_task(void *param)
{
#if (ENABLE_MAGNETOMETER != 0)
	ist8130_init();
	freertos_task_delay(10);
#endif

#if (ENABLE_RANGEFINDER != 0)
	lidar_lite_init();
	freertos_task_delay(10);
#endif

	while(ins_sync_buffer_is_ready() == false);

	while(1) {
#if ((ENABLE_MAGNETOMETER != 0) && (ENABLE_RANGEFINDER != 0))
		ist8310_read_sensor();
		lidar_lite_read_sensor();
		lidar_lite_read_sensor();
#elif (ENABLE_MAGNETOMETER != 0)
		ist8310_read_sensor();
#elif (ENABLE_RANGEFINDER != 0)
		lidar_lite_read_sensor();
#endif
	}
}

void f4_sw_i2c_driver_register_task(const char *task_name, configSTACK_DEPTH_TYPE stack_size,
                                    UBaseType_t priority)
{
	xTaskCreate(f4_sw_i2c_driver_task, task_name, stack_size, NULL, priority, NULL);
}

