const int inPin = A9;
int a;
char ch;
long int t0, t1;
int state;
String line = "";

// the setup() method runs once, when the sketch starts

void setup() {
  // initialize the digital pin as an output.
  Serial.begin(9600);
  pinMode(inPin, INPUT);
  state = 0;
}

// the loop() methor runs over and over again,
// as long as the board has power

void loop() {
  if(Serial.available()){
    ch = Serial.read();
    if (ch == '\n') {
      //Serial.println(line);
      if (line.equals("read")){
        a = analogRead(inPin);   // set the LED on
        Serial.println(a);
      }
      line = "";
    } else if (ch != '\r') {
      line += ch;  
    }
  }
}
