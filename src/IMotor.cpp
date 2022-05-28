
#include "IMotor.h"

IMotor::IMotor(uint8_t IN1, uint8_t IN2, uint8_t pwmpin, int lower_limit, int upper_limit)
{
  this->IN1 = IN1;
  this->IN2 = IN2;
  _pwmpin = pwmpin;
  _upper_limit = upper_limit;
  _lower_limit = lower_limit;
}

// IMotor::IMotor(uint8_t IN, uint8_t pwmpin, int lower_limit, int upper_limit)
// {
//   this->IN1 = IN;
//   this->IN2 = -1;
//   _pwmpin = pwmpin;
//   _upper_limit = upper_limit;
//   _lower_limit = lower_limit;
// }

void IMotor::init(int _pwmChannel)
{
  pinMode(_pwmpin, OUTPUT);
  pinMode(IN1, OUTPUT);
  pinMode(IN2, OUTPUT);
  motor_state = 1;
  pwmChannel = _pwmChannel;
  // configure LED PWM functionalitites
  ledcSetup(pwmChannel, freq, resolution);
  // attach the channel to the GPIO to be controlled
  ledcAttachPin(_pwmpin, pwmChannel);
}

void IMotor::setMotor(int pwmVal)
{
  setMotor(sgn(pwmVal), abs(pwmVal));
}

void IMotor::setMotor(int dir, int pwmVal)
{
  pwmVal = constrain(pwmVal, _lower_limit, _upper_limit);

  if (dir == 1 && motor_state == 1)
  {
    digitalWrite(IN1, HIGH);
    digitalWrite(IN2, LOW);
  }
  else if (dir == -1 && motor_state == 1)
  {
    digitalWrite(IN1, LOW);
    digitalWrite(IN2, HIGH);
  }
  else
  {
    digitalWrite(IN1, LOW);
    digitalWrite(IN2, LOW);
  }

  if (_pwmpin != 0)
  {
    ledcWrite(pwmChannel, pwmVal);
  }
}

void IMotor::setMotor(int dir, int pwmVal, int in1, int in2)
{
  pwmVal = constrain(pwmVal, _lower_limit, _upper_limit);
  if (_pwmpin != 0)
    ledcWrite(pwmChannel, pwmVal);
 
  if (dir == 1 && motor_state == 1)
  {
    digitalWrite(in1, HIGH);
    digitalWrite(in2, LOW);
  }
  else if (dir == -1 && motor_state == 1)
  {
    digitalWrite(in1, LOW);
    digitalWrite(in2, HIGH);
  }
  else
  {
    digitalWrite(in1, LOW);
    digitalWrite(in2, LOW);
  }
}

void IMotor::turn_on()
{
  motor_state = 1;
}

void IMotor::turn_off()
{
  motor_state = 0;
  if (pwmChannel != 0)
    ledcWrite(pwmChannel, 0);
  (_pwmpin != 0) ? digitalWrite(IN1, 0) : ledcWrite(IN1, 0);
  (_pwmpin != 0) ? digitalWrite(IN2, 0) : ledcWrite(IN2, 0);
}

void IMotor::limit(int lower_limit, int upper_limit)
{
  _lower_limit = lower_limit;
  _upper_limit = upper_limit;
}

bool IMotor::target_reached(bool reset)
{
  if (reset)
    target_is_reached = false;
  return target_is_reached;
}
