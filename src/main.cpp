/*
 * Blink
 * Turns on an LED on for one second,
 * then off for one second, repeatedly.
 */

#include <Arduino.h>
#include "ESP32Encoder.h"
#include <InterruptEncoder.h>
#include <WiFi.h>
#include <WiFiUdp.h>

#include "IMotor.h"
#include "LQRMotors.h"
#include "eeprom_utils.h"

#include "IMath/IMaths.h"

#define pushButton_pin 18 
#define PIN_SDA 21
#define PIN_SDL 22

MPU9250 mpu;

// Структора данных , кнтроля робокара
struct dataControl
{
  int speedX; // Скорость вращения мотрчиков справа , ШИМ-А
  int speedY; // Скорость вращения мотрчиков слева  , ШИМ-B
  int turn;
  bool is_null_pos;

  float kp;
  float kd;
  float ki;
  float kf;
  float kt;
};

static dataControl *data_control =  nullptr;
static IMath::Vector3 euler;

/**/
struct Quaternion
{
  float w, x, y, z;
};

struct DataDescriptor
{
  int num;
  Quaternion Quat;
  float EulerX;
  float EulerY;
  float EulerZ;

  //--------------//

  float AccBiasX;
  float AccBiasY;
  float AccBiasZ;

  float GyroBiasX;
  float GyroBiasY;
  float GyroBiasZ;

  float MagBiasX;
  float MagBiasY;
  float MagBiasZ;

  float MagScaleX;
  float MagScaleY;
  float MagScaleZ;

  //--------------//

} _dataDescriptor;
/**/

// void print_roll_pitch_yaw()
// {
//     Serial.print("Yaw, Pitch, Roll: ");
//     Serial.print(mpu.getYaw(), 2);
//     Serial.print(", ");
//     Serial.print(mpu.getPitch(), 2);
//     Serial.print(", ");
//     Serial.println(mpu.getRoll(), 2);
// }

void print_roll_pitch_yaw();

void print_calibration()
{
  Serial.println("< calibration parameters >");
  Serial.println("accel bias [g]: ");
  Serial.print(mpu.getAccBiasX() * 1000.f / (float)MPU9250::CALIB_ACCEL_SENSITIVITY);
  Serial.print(", ");
  Serial.print(mpu.getAccBiasY() * 1000.f / (float)MPU9250::CALIB_ACCEL_SENSITIVITY);
  Serial.print(", ");
  Serial.print(mpu.getAccBiasZ() * 1000.f / (float)MPU9250::CALIB_ACCEL_SENSITIVITY);
  Serial.println();
  Serial.println("gyro bias [deg/s]: ");
  Serial.print(mpu.getGyroBiasX() / (float)MPU9250::CALIB_GYRO_SENSITIVITY);
  Serial.print(", ");
  Serial.print(mpu.getGyroBiasY() / (float)MPU9250::CALIB_GYRO_SENSITIVITY);
  Serial.print(", ");
  Serial.print(mpu.getGyroBiasZ() / (float)MPU9250::CALIB_GYRO_SENSITIVITY);
  Serial.println();
  Serial.println("mag bias [mG]: ");
  Serial.print(mpu.getMagBiasX());
  Serial.print(", ");
  Serial.print(mpu.getMagBiasY());
  Serial.print(", ");
  Serial.print(mpu.getMagBiasZ());
  Serial.println();
  Serial.println("mag scale []: ");
  Serial.print(mpu.getMagScaleX());
  Serial.print(", ");
  Serial.print(mpu.getMagScaleY());
  Serial.print(", ");
  Serial.print(mpu.getMagScaleZ());
  Serial.println();
}

//---------------------------------//

unsigned int localPort = 8888; // Локальный порт прослушевания сети

/**
const char *ssid = "wera";           // SSID имя WiFi точки-доступа робокара
const char *password = "123qwe1023"; // Пароль WiFii

IPAddress local_ip(192,168,1,100); // IP-адрес робокара
IPAddress gateway(192,168,1,1);    // IP-адрес шлюза  
IPAddress subnet(255,255,255,0);   // Подсеть 
/**/

/**/
const char *ssid     = "ROBO_CAR"; // SSID имя WiFi точки-доступа робокара 
const char* password = "12345678"; // Пароль WiFii  

IPAddress local_ip(192,168,43,100); // IP-адрес робокара
IPAddress gateway(192,168,1,1);    // IP-адрес шлюза  
IPAddress subnet(255,255,255,0);    // Подсеть 
/**/

#define PACKET_MAX_SIZE 255         // Масимальный размер пакета-данных
char packetBuffer[PACKET_MAX_SIZE]; // Буферы для приема и отправки пакета-данных,
char replayBuffer[PACKET_MAX_SIZE];

WiFiUDP Udp; // UPD-сокет обект
//----------------------------------//



//----------------------------------//

#define PIN_IN1 19
#define PIN_IN2 18 

#define PIN_IN3 35
#define PIN_IN4 34

// Setup a RotaryEncoder with 2 steps per latch for the 2 signal input pins:
InterruptEncoder encoderA;
InterruptEncoder encoderB;

// //--------------------------//

// Motor A
int motor1Pin1 = 27;
int motor1Pin2 = 26;
int enable1Pin = 14;

// Motor B
int motor2Pin1 = 25;
int motor2Pin2 = 33;
int enable2Pin = 32;

#define LEFT_PWM_MIN 0//50//35  
#define RIGHT_PWM_MIN 0//50//35 
#define LEFT_MAX_SPEED_FORCE 250//200//55
#define RIGHT_MAX_SPEED_FORCE 250//200//55

IMotor motorA(motor1Pin1, motor1Pin2, enable1Pin, LEFT_PWM_MIN, LEFT_MAX_SPEED_FORCE);
IMotor motorB(motor2Pin1, motor2Pin2, enable2Pin, RIGHT_PWM_MIN, RIGHT_MAX_SPEED_FORCE);

int dutyCycle = 180;

// float k[4] ={0.2 ,  -7746.3  , 0.025 ,  -70};//22 ,-7458.4,  40  ,-60.8};
float k[4] = {18.8, -7746.3, 40.1, -70}; // 22 ,-7458.4,  40  ,-60.8};

//  gains for yaw | yaw_dot
// float k1[2] = {5000, 500};
float k1[2] = {4000, 80.03};

LQRMotors LQRRobot(&motorA, &motorB, &encoderA, &encoderB);

void setup()
{
  Serial.begin(115200);
  //-----------------------------------------------//
  Wire.begin(PIN_SDA, PIN_SDL);
  // Wire.setClock(400000);
  delay(1000);

  MPU9250Setting setting;
  setting.accel_fs_sel = ACCEL_FS_SEL::A16G;
  setting.gyro_fs_sel = GYRO_FS_SEL::G2000DPS;
  setting.mag_output_bits = MAG_OUTPUT_BITS::M16BITS;
  setting.fifo_sample_rate = FIFO_SAMPLE_RATE::SMPL_200HZ;
  setting.gyro_fchoice = 0x03;
  setting.gyro_dlpf_cfg = GYRO_DLPF_CFG::DLPF_41HZ;
  setting.accel_fchoice = 0x01;
  setting.accel_dlpf_cfg = ACCEL_DLPF_CFG::DLPF_45HZ;

  if (!mpu.setup(0x68, setting))
  {
    // change to your own address
    while (true)
    {
      Serial.println("MPU connection failed. Please check your connection with `connection_check` example.");
      delay(5000);
    }
  }

#if defined(ESP_PLATFORM) || defined(ESP8266)
  EEPROM.begin(0x80);
#endif

  // load from eeprom
  loadCalibration();

  // mpu.setFilterIterations(5);
  // mpu.selectFilter(QuatFilterSel::MADGWICK);

  mpu.setFilterIterations(10);
  mpu.selectFilter(QuatFilterSel::MAHONY);

  ////mpu.Init();

  delay(3000);

  //---------------------//

  encoderA.attach(PIN_IN1, PIN_IN2);
  encoderB.attach(PIN_IN3, PIN_IN4);

  //---------------------//

  motorA.init(1);
  motorB.init(0);

  //---------------------//

  LQRRobot.Init();

  //--------------------//

  // присваиваем статичесий IP адрес
  WiFi.mode(WIFI_STA); // режим клиента
  WiFi.config(local_ip, gateway, subnet);
  WiFi.begin(ssid, password);

  int n = 0;
  while (WiFi.status() != WL_CONNECTED)
  {
    Serial.print('.');
    delay(500);

    if (n++ > 20)
      break;
  }

  if (WiFi.status() != WL_CONNECTED)
  {
    boolean result =
    WiFi.softAP(ssid, password);                  // Устанавливаем режым точки доступа WiFi
    WiFi.softAPConfig(local_ip, gateway, subnet); // Устанавливаем статические IP-адреса
    delay(100);                                   // Ждем 100 милисекунд

    // вывод данных о сети
    IPAddress myIP = WiFi.softAPIP(); // IP-адрес робокара
    Serial.print("AP IP address: ");
    Serial.println(myIP);
    if (result == true)
      Serial.println("Ready");
    else
      Serial.println("Failed!");
  }

  // выввод информации о сервере
  Serial.print("Connected! IP address: ");
  Serial.println(WiFi.localIP());
  Serial.printf("UDP server on port %d\n", localPort);

  Udp.begin(localPort); // Начинаем слушать порт 8888 . Ждем потключения клиента

  LQRRobot.AbsolutResetEncoderPos();
  
  LQRRobot.pwm_right_offset = 0;
  LQRRobot.pwm_left_offset = 0;

  LQRRobot.position_set_point = 0;
  LQRRobot.velocity_set_point = 0;

  LQRRobot.position_set_point = 0;
  LQRRobot.average_theta = 0;

}

// Read the current position of the encoder and print out when changed.
void loop()
{
  if (mpu.update())
  {
    static uint32_t prev_ms = millis();
    if (millis() > prev_ms + 25)
    {
      print_roll_pitch_yaw();
      prev_ms = millis();
    }

    /**

    // Ждем ..! И если есть данные, начинаем обрабатывать пакет-данных
    int packetSize = Udp.parsePacket();
    if (packetSize)
    {
      // Читаем пакеты в packageBufffer
      int n = Udp.read(packetBuffer, sizeof(dataControl));
      packetBuffer[n] = 0;

      // Переобразуем данные в удобный нам формат
      dataControl *data_control = (dataControl *)&packetBuffer;

      // data_control->speedX = map(data_control->speedX,0,1024,0,255);
      // data_control->speedY = map(data_control->speedY,0,1024,0,255);

      // if (data_control->is_null_pos)
      // {
      //   LQRRobot.AbsolutResetEncoderPos();
      // }

      k[0] = data_control->kp;
      k[2] = data_control->kd;

      k1[0] = data_control->ki;
      k1[1] = data_control->kf;

      LQRRobot.position_set_point += data_control->speedX / 100.f;
      LQRRobot.pwm_right_offset = -data_control->turn / 20.f;
      LQRRobot.pwm_left_offset = data_control->turn / 20.f;


      LQRRobot.offset_Yaw -= data_control->turn / 500.f;


      IMath::Quaternion Q(mpu.getQuaternionW(),
                          mpu.getQuaternionX(),
                          mpu.getQuaternionY(),
                          mpu.getQuaternionZ());

      Q = Q * IMath::Quaternion::FromAxisRot(IMath::Vector3::Z , LQRRobot.offset_Yaw);
      Q.Normalize();

      // if(abs(pwm_left_offset) > 0 && abs(pwm_left_offset) < 50) pwm_left_offset = 50 * sgn(pwm_left_offset);
      // if(abs(pwm_right_offset) > 0 && abs(pwm_right_offset) < 50) pwm_right_offset = 50 * sgn(pwm_right_offset);

      Serial.println(Q.GetEulerAngles3().x);
      Serial.println(Q.GetEulerAngles3().y);
      Serial.println(Q.GetEulerAngles3().z);

      //Serial.println(mpu.getYaw());
      LQRRobot.update_encoder_states(Q.GetEulerAngles3().z);
      LQRRobot.update_motors(data_control->is_null_pos, data_control->speedX, k, k1, data_control->ki);
      LQRRobot.SerialPlotterPrint();

      //  char buff[256];
      // sprintf( buff , "sx:%d , sy:%d  \n" ,
      //          data_control->speedX ,
      //          data_control->speedY );

      char buff[256];
      sprintf(buff, "sx:%d , sy:%d  \n",
              encoderA.read(),
              encoderB.read());

      // Выввод присланых данных в сериал порт
      // Serial.printf(buff);

      Udp.beginPacket(Udp.remoteIP(), Udp.remotePort());
      Udp.write((uint8_t *)&buff, 255);
      Udp.endPacket();

     // delay(50);
    }


    **/

    // k[0] = 10.0;//data_control->kp;
    // k[2] = 0.2;//data_control->kd;

    // LQRRobot.position_set_point = 1000;
    // // LQRRobot.pwm_right_offset = -data_control->turn / 25;
    // // LQRRobot.pwm_left_offset = data_control->turn / 25;

    // // if(abs(pwm_left_offset) > 0 && abs(pwm_left_offset) < 50) pwm_left_offset = 50 * sgn(pwm_left_offset);
    // // if(abs(pwm_right_offset) > 0 && abs(pwm_right_offset) < 50) pwm_right_offset = 50 * sgn(pwm_right_offset);

    // LQRRobot.update_encoder_states();
    // LQRRobot.update_motors(false, 0, k, k1, 0);
    // LQRRobot.SerialPlotterPrint();
  }

  // delay(50);
}

void print_roll_pitch_yaw()
{
  // Ждем ..! И если есть данные, начинаем обрабатывать пакет-данных
  int packetSize = Udp.parsePacket();
  if (packetSize)
  {
    // Читаем пакеты в packageBufffer
    int n = Udp.read(packetBuffer, sizeof(dataControl));
    packetBuffer[n] = 0;

    // Переобразуем данные в удобный нам формат
    data_control = (dataControl *)&packetBuffer;

    // data_control->speedX = map(data_control->speedX,0,1024,0,255);
    // data_control->speedY = map(data_control->speedY,0,1024,0,255);

    // if (data_control->is_null_pos)
    // {
    //   LQRRobot.AbsolutResetEncoderPos();
    // }

    k[0] = data_control->kp;
    k[2] = data_control->kd;

    k1[0] = data_control->ki;
    k1[1] = data_control->kf;



    LQRRobot.position_set_point += data_control->speedX / 100.f;
    LQRRobot.pwm_right_offset = -data_control->turn / 20.f;
    LQRRobot.pwm_left_offset = data_control->turn / 20.f;

    if (IMath::abs(data_control->turn) > 0.0001f)
    {
      LQRRobot.offset_Yaw -= data_control->turn / 500.f;
    }
      

    IMath::Quaternion Q(mpu.getQuaternionW(),
                        mpu.getQuaternionX(),
                        mpu.getQuaternionY(),
                        mpu.getQuaternionZ());

    IMath::Quaternion Q_Target = 
    IMath::Quaternion::FromAxisRot(IMath::Vector3::Z, LQRRobot.offset_Yaw);
    Q = Q * Q_Target;
    Q.Normalize();

    // if(abs(pwm_left_offset) > 0 && abs(pwm_left_offset) < 50) pwm_left_offset = 50 * sgn(pwm_left_offset);
    // if(abs(pwm_right_offset) > 0 && abs(pwm_right_offset) < 50) pwm_right_offset = 50 * sgn(pwm_right_offset);

    // Serial.println(Q.GetEulerAngles3().x);
    // Serial.println(Q.GetEulerAngles3().y);
    // Serial.println(Q.GetEulerAngles3().z);

    euler = Q.GetEulerAngles3();
    IMath::Vector3 euler_target = Q_Target.GetEulerAngles3();


    // Serial.println("Real,Target");
    // Serial.print(/*IMath::IRadiansToDegrees*/(mpu.getEulerZ()));
    // Serial.print("  ");
    // Serial.print(euler_target.z);
    // Serial.println();


    Serial.println("Control_Position,Update_Position");
    Serial.print(/*IMath::IRadiansToDegrees*/(LQRRobot.position_set_point));
    Serial.print("  ");
    Serial.print(LQRRobot.average_theta);
    Serial.println();
    

    // Serial.print("Yaw, Pitch, Roll: ");
    // Serial.print(euler.z, 2);
    // Serial.print(", ");
    // Serial.print(mpu.getPitch(), 2);
    // Serial.print(", ");
    // Serial.println(mpu.getRoll(), 2);

    //  char buff[256];
    // sprintf( buff , "sx:%d , sy:%d  \n" ,
    //          data_control->speedX ,
    //          data_control->speedY );

    /**/


    _dataDescriptor.num = 20;
    _dataDescriptor.Quat.x = mpu.getQuaternionX();
    _dataDescriptor.Quat.y = mpu.getQuaternionY();
    _dataDescriptor.Quat.z = mpu.getQuaternionZ();
    _dataDescriptor.Quat.w = mpu.getQuaternionW();

    _dataDescriptor.EulerX = mpu.getEulerX();
    _dataDescriptor.EulerY = mpu.getEulerY();
    _dataDescriptor.EulerZ = mpu.getEulerZ();

    char buff[256];
    sprintf(buff, "sx:%d , sy:%d  \n",
            encoderA.read(),
            encoderB.read());

    // // Выввод присланых данных в сериал порт
    // // Serial.printf(buff);

   if(data_control)
   {
        // Serial.println(mpu.getYaw());
      LQRRobot.update_encoder_states(euler.z);
      LQRRobot.update_motors(data_control->is_null_pos,
                             data_control->speedX, k, k1, 
                             data_control->ki);

                             data_control->turn = data_control->speedX = 0;
      //LQRRobot.SerialPlotterPrint();
   }

    Udp.beginPacket(Udp.remoteIP(), Udp.remotePort());
    Udp.write((uint8_t*)&_dataDescriptor, sizeof(DataDescriptor));
    Udp.endPacket();

    // delay(50);
  }
  else
  {
    /**
    LQRRobot.AbsolutResetEncoderPos();
    LQRRobot.position_set_point = 0;
    LQRRobot.average_theta = 0;
    /**/
  }


  //  if(data_control)
  //  {
  //       // Serial.println(mpu.getYaw());
  //     LQRRobot.update_encoder_states(euler.z);
  //     LQRRobot.update_motors(data_control->is_null_pos,
  //                            data_control->speedX, k, k1, 
  //                            data_control->ki);

  //                            data_control->speedX=0;
  //     //LQRRobot.SerialPlotterPrint();
  //  }
   
}