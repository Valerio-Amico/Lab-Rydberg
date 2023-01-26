'''
This script sends a periodic serial request to arduino, for measuring the a tension,
and saves all the data in a file .txt.
'''

import serial
import time
from os import path

def now(str):
    return time.strftime(str,time.gmtime(time.time()))

serial_port = '/dev/tty.usbmodem123383401'
baud_rate = 9600 #In arduino, Serial.begin(baud_rate)

delay = 0.1 # seconds
today = time.strftime("%d-%m-%Y",time.gmtime(time.time()))
name_exp = "prova_ripompa_before_and_after_with_lamina_good"
name = "dati/output"+today+name_exp+".txt"

if path.exists(name)==True:
    output_file = open(name, "a")
else:
    output_file = open(name, "w")

ser = serial.Serial(serial_port, baud_rate, timeout=0.1)

while True:
    # write with the serial connection to arduino
    if ser.is_open==False:
        ser.open()
    ser.write("read\n".encode('utf-8'))
    # read the temperature from arduino
    # time.sleep(0.05)
    while ser.in_waiting==False:
        pass
    line = ser.readline()
    line = line.decode("utf-8") #ser.readline returns a binary, convert to string
    #else:
    #    line = "nothing"

    if ser.is_open==True:
        ser.close()
    
    print(line)
    
    # Getting current date and time
    output_file = open(name, "a")

    output_file.write(now("%d-%m-%Y %H:%M:%S"))
    output_file.write(" ")
    output_file.write(line)

    output_file.close()
    print(now("%d-%m-%Y %H:%M:%S"), " ", line)

    time.sleep(delay)

