##############
## Script listens to serial port and writes contents into a file
##############
## requires pySerial to be installed 

import serial
import datetime
from os import path
  
serial_port = '/dev/tty.usbmodem1101'
baud_rate = 9600 #In arduino, Serial.begin(baud_rate)

today = datetime.datetime.now()
name = "arduino_for_conditioner/dati/output"+today.strftime("%m-%d-%Y")+".txt"

print(path.exists(name))

if path.exists(name)==True:
    output_file = open(name, "a")
else:
    output_file = open(name, "w")

ser = serial.Serial(serial_port, baud_rate)

while True:
    line = ser.readline()
    line = line.decode("utf-8") #ser.readline returns a binary, convert to string
    
    # Getting current date and time
    now = datetime.datetime.now()
    output_file.write(now.strftime("%m-%d-%Y %H:%M:%S"))
    output_file.write(" ")
    output_file.write(line)
    print(now.strftime("%m-%d-%Y"), " ", line, "\n")