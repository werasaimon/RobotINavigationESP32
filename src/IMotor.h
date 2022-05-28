#ifndef MOTOR_H
#define MOTOR_H


#include <Arduino.h>

/*--------------------------------------------*/
class IMotor
{
  
public:

  bool motor_state; //this variable let the pwm signal go
  IMotor(uint8_t IN1, uint8_t IN2, uint8_t pwmpin = 0, int lower_limit = 50, int upper_limit = 255); // constructor
  //IMotor(uint8_t IN, uint8_t pwmpin = 0, int lower_limit = 50, int upper_limit = 255); // constructor
  ~IMotor(){};

  void setMotor(int pwmVal);
  void setMotor(int dir, int pwmVal);
  

  void init(int _pwmChannel); // this function initializes pid regulator
  void turn_on();  // changes motor_state variable to true
  void turn_off(); // changes motor_state variable to false
  

  void limit(int,int);
  bool target_reached(bool reset=false);// check if target position is reached by motor

  //-----------------------------//

  uint8_t IN2, IN1;
  void setMotor(int dir, int pwmVal, int in1, int in2);


  void setRange(int upper_limit, int lower_limit )
  {
    _upper_limit = upper_limit;
    _lower_limit = lower_limit;
  } 
 
 public:
  
  bool buffer;
  uint8_t _pwmpin = 0;
  int _upper_limit = 0, _lower_limit = 0;
  bool target_is_reached=false;

  template <typename T> int sgn(T val) 
  {
      return (T(0) < val) - (val < T(0));
  }

  int pwmChannel = 0;
  const int freq = 30000;
  const int resolution = 8;

};


#endif // MOTOR_H