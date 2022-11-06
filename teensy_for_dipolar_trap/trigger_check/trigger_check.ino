const int outPin = 13;
const int inPin = 23;
int a;
long int t0, t1;
int state;

// the setup() method runs once, when the sketch starts

void setup() {
  // initialize the digital pin as an output.
  Serial.begin(9600);
  pinMode(outPin, OUTPUT);
  pinMode(inPin, INPUT);
  state = 0;
}

// the loop() methor runs over and over again,
// as long as the board has power

void loop() {
  //t0 = micros(); 
  a = digitalRead(inPin);   // set the LED on
  if(a==1){
    digitalWrite(outPin, HIGH);
    //state=1;
    //t1 = micros();
    //t1 = t1-t0;
    //Serial.println(a);
    //Serial.print("tempo [us]: ");
    //Serial.println(t1);
    delayNanoseconds(5000);
    digitalWrite(13, LOW);
    delayMicroseconds(10);
    //state=0;
  }
  //if(a==0 and state==1){
  //  digitalWrite(13, LOW);
  //  state=0;
  //}
}
