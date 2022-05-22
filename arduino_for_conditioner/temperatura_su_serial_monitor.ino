#include <OneWire.h> 
#include <DallasTemperature.h>

#define ONE_WIRE_BUS 10  // digital pin for Dallas temperature sensor

OneWire oneWire(ONE_WIRE_BUS);
DallasTemperature sensors(&oneWire);

double Temperature;
char ch;
String line = "";
bool FIRST = true;

void setup() {
  // initialize serial communication at 9600 bits per second:
  Serial.begin(9600);
  sensors.begin();
}

// the loop routine runs over and over again forever:
void loop() {

  if(Serial.available()){
    ch = Serial.read();
    if (ch == '\n') {
      Serial.println(line);
      if (line.equals("read")){
        sensors.requestTemperatures();
        Temperature=sensors.getTempCByIndex(0);
        Serial.println(Temperature);
      }
      line = "";
    } else if (ch != '\r') {
      line += ch;  
    }
  }
}