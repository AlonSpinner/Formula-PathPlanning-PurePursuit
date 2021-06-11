import serial
ArduinoData=serial.Serial('COM3',9600)

SteeringAngle=24
ArduinoData.write(chr(SteeringAngle))