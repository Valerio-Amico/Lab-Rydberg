#include <OneWire.h> 
#include <DallasTemperature.h>

#define ONE_WIRE_BUS 10  // digital pin for Dallas temperature sensor

OneWire oneWire(ONE_WIRE_BUS);
DallasTemperature sensors(&oneWire);

double Temperature;

long unsigned delay_time = 60*1000;
long unsigned t0 = millis();

void setup() {
  // initialize serial communication at 9600 bits per second:
  Serial.begin(9600);
  sensors.begin();
}

// the loop routine runs over and over again forever:
void loop() {

  sensors.requestTemperatures();
  Temperature=sensors.getTempCByIndex(0);
  Serial.println(Temperature);

  while ((millis() - t0) < delay_time){
    delay(1);
  }
  t0 = millis();
}