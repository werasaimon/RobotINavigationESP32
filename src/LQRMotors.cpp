
#include "LQRMotors.h"
// #include "eeprom_utils.h"


//MPU9250 mpu;

#define MIN_EPS 0.01

LQRMotors::LQRMotors(IMotor *m1, IMotor *m2, InterruptEncoder *en1, InterruptEncoder *en2)
{
    m_Motors[0] = m1;
    m_Motors[1] = m2;
    m_Encoders[0] = en1;
    m_Encoders[1] = en2;
}

LQRMotors::~LQRMotors()
{
}

void LQRMotors::Init()
{
    offset_Yaw = 0;
    position_set_point = 0;
    velocity_set_point = 0;
    left_encoder_count = 0;
    right_encoder_count = 0;
    left_RPM = 0;
    right_RPM = 0;
    left_prev_count = 0;
    right_prev_count = 0;

    average_theta = 0;
    average_RPM = 0;

    v[2] = {0};
    u[2] = {0};

    yaw = 0;
    previous_yaw = 0;
    yaw_dot = 0;

    m_Motors[0]->init(0);
    m_Motors[1]->init(1);

    //------------------//

//     if (!mpu.setup(0x68))
//     { // change to your own address
//         while (true)
//         {
//             Serial.println("MPU connection failed. Please check your connection with `connection_check` example.");
//             delay(5000);
//         }
//     }

// #if defined(ESP_PLATFORM) || defined(ESP8266)
//     EEPROM.begin(0x80);
// #endif

//     // load from eeprom
//     loadCalibration();

//     mpu.setFilterIterations(4);
//     mpu.selectFilter(QuatFilterSel::MADGWICK);

    //------------------//
}

void LQRMotors::SetDispatch(const LQR_Dispatch &_Dispatch)
{
    m_Dispatch = _Dispatch;
}

void LQRMotors::update_encoder_states(float mpuYaw)
{
    //----------------------------------------------------//

    // Make a local copy of the global encoder count
    volatile float left_current_count = m_Encoders[0]->read();
    volatile float right_current_count = m_Encoders[1]->read();

    //                (Change in encoder count)
    // RPS =   __________________________________________
    //          (Change in time --> 5ms) * (PPR --> 540)
    left_RPM = (float)((left_current_count - left_prev_count) * sRPM_radian_converter);
    right_RPM = (float)((right_current_count - right_prev_count) * sRPM_radian_converter);

    // Store current encoder count for next iteration
    left_prev_count = left_current_count;
    right_prev_count = right_current_count;

    // previous_yaw = yaw;
    // yaw = (+right_current_count - left_current_count) * sYaw_radian_multiplier;
    // yaw_dot = (yaw - previous_yaw) * sone_by_dT;

    average_theta = ((left_current_count + right_current_count) * 0.0116f); // 2pi/540
    average_RPM = (left_RPM + right_RPM) * 0.5f;

    //----------------------------------------------------//

    previous_yaw = yaw;
    yaw = mpuYaw * MIN_EPS;// * sYaw_radian_multiplier;
    yaw_dot = (yaw - previous_yaw) * 100;//MIN_EPS;// * sone_by_dT;

    //----------------------------------------------------//
}

void LQRMotors::update_motors(bool isNull, int vel, float k[4], float k1[2], int min_dist)
{
    // error
    float e = (position_set_point - average_theta);

    if (fabs(e) < 20.f)
    {
      e = 0;
    }
       
    error = e;

    u[0] = (k[0] * e + k[2] * (velocity_set_point - average_RPM));
    u[1] = (k1[0] * yaw + k1[1] * yaw_dot);
    v[0] = (0.5f * (u[0] + u[1])); // motor left
    v[1] = (0.5f * (u[0] - u[1])); // motor right

    // v[0] = (0.5f * u[0]); // motor left
    // v[1] = (0.5f * u[0]); // motor right


    if (/*fabs(e) < min_dist && */isNull)
    {
        // v[0]=v[1]=0;
        //if (vel == 0)
        {
            m_Encoders[0]->count = ((position_set_point / 0.0116f) * -0.5f);
            m_Encoders[1]->count = ((position_set_point / 0.0116f) * -0.5f);

            velocity_set_point = average_RPM = 0;
            position_set_point = average_theta = 0;
        }
    }

    v[0] = v[0] + pwm_left_offset;
    v[1] = v[1] + pwm_right_offset;
    // v[0]=(fabs(v[0]) < LEFT_PWM_MIN) ? LEFT_PWM_MIN * sgn(v[0]): v[0];
    // v[1]=(fabs(v[1]) < RIGHT_PWM_MIN) ? RIGHT_PWM_MIN * sgn(v[1]): v[1];

    m_Motors[0]->setMotor(v[0]);
    m_Motors[1]->setMotor(v[1]);


       // Serial.println(error);
    // Serial.println("Yaw , Yaw_Dot");
    // Serial.print(k1[0] * yaw);
    // Serial.print("  ");
    // Serial.print(k1[1] * yaw_dot);
    // Serial.println();
}

void LQRMotors::SerialPlotterPrint()
{
    /**
    Serial.println("position_set_point,average_theta,error,Yaw");
    Serial.print(position_set_point);
    Serial.print("  ");
    Serial.print(average_theta);
    Serial.print("  ");
    Serial.print(error);
    Serial.print("  ");
    Serial.print(yaw);
    Serial.println();
    /**/

    Serial.println("Yaw,Yaw_Dot");
    Serial.print(yaw);
    Serial.print("  ");
    Serial.print(yaw_dot);
    Serial.println();
}

void LQRMotors::AbsolutResetEncoderPos()
{
    m_Encoders[0]->count = ((position_set_point / 0.0116) * -0.5);
    m_Encoders[1]->count = ((position_set_point / 0.0116) * -0.5);
}