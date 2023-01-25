const int inPin1 = A9;
const int inPin2 = A8;
int a, b;
char ch;
long int t0, t1;
int state;
String line = "";

// the setup() method runs once, when the sketch starts

void setup() {
  // initialize the digital pin as an output.
  Serial.begin(9600);
  pinMode(inPin1, INPUT);
  pinMode(inPin2, INPUT);
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
        a = analogRead(inPin1); 
        b = analogRead(inPin2);
        Serial.print(a);
        Serial.print(" ");
        Serial.println(b);
      }
      line = "";
    } else if (ch != '\r') {
      line += ch;  
    }
  }
}
