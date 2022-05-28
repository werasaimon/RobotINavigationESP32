#ifndef LQRMOTORS_H
#define LQRMOTORS_H

#include "InterruptEncoder.h"
#include <MPU9250.h>
#include "IMotor.h"


#define sRPM_radian_converter 13.962//2.327
#define sYaw_radian_multiplier 0.0020945
#define sencoder_increment_rate 3

//float dT   = 0.02;              //0.012
#define sdT 0.005
#define sone_by_dT 200
#define salpha1 0.03

// #define LEFT_PWM_MIN 50//60//45
// #define RIGHT_PWM_MIN 50//50//70

struct LQR_Dispatch
{
    volatile float left_encoder_count;
    volatile float right_encoder_count;
    float encoder_set_point;
    float velocity_set_point;

    float pwm_right_offset;
    float pwm_left_offset;
    float average_theta;
    float average_RPM;
    float v[2];
    float u[2];

    float yaw;
    float previous_yaw;
    float yaw_dot;
};

class LQRMotors
{

public:
 
    float offset_Yaw;

    float pwm_right_offset = 0;
    float pwm_left_offset = 0;

    float position_set_point = 0;
    float velocity_set_point = 0;
//private:
   public: 
    /* data */


    volatile float error;
    volatile float left_encoder_count = 0;
    volatile float right_encoder_count = 0;
    volatile float left_RPM = 0;
    volatile float right_RPM = 0;
    volatile float left_prev_count = 0;
    volatile float right_prev_count = 0;

    float average_theta = 0;
    float average_RPM = 0;

    float v[2] = {0};
    float u[2] = {0};

    float yaw = 0;
    float previous_yaw = 0;
    float yaw_dot = 0;

    LQR_Dispatch m_Dispatch;

    InterruptEncoder *m_Encoders[2];
    IMotor *m_Motors[2];

public:
     LQRMotors(IMotor *m1, IMotor *m2, InterruptEncoder *en1, InterruptEncoder *en2);
    ~LQRMotors();

    void Init();
    void SetDispatch(const LQR_Dispatch &_Dispatch);

    void update_encoder_states(float mpuYaw);
    void update_motors(bool isNull, int vel,float k[4],float k1[2], int min_dist);

    void AbsolutResetEncoderPos();

    void SerialPlotterPrint();

};

#endif // LQRMOTORS_H