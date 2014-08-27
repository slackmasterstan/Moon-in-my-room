#include <stdlib.h>
#include <math.h>
#include <Wire.h>
#include "RTClib.h"


// South and West are negative ! North and East are positive !
#define LONGITUDE -97.2672
#define LATITUDE  35.4825
#define TZOFFSET -5

#define LED1 4
#define LED2 9
#define LED3 10
#define LED4 12
#define LED5 14
#define LED6 15

#define LEDMAX 255

RTC_DS1307 RTC;

static float Dec[3] = { 0.0, 0.0, 0.0 };
static float RAn[3] = { 0.0, 0.0, 0.0 };
static float VHz[3] = { 0.0, 0.0, 0.0 };
static float Sky[3] = { 0.0, 0.0, 0.0 };
const static float DR = M_PI / 180;

static int Rise_time[2] = {0.0, 0.0}; 
static int Set_time[2] = {0.0, 0.0};

static bool Moonrise, Moonset;
int then;
int start=1;

// returns value for sign of argument
float sgn(float x)
{
    float rv;
    if (x > 0.0)
	rv = 1;
    else if (x < 0.0)
	rv = -1;
    else
	rv = 0;
    return rv;
}

// determine Julian day from calendar date
// (Jean Meeus, "Astronomical Algorithms", Willmann-Bell, 1991)
float julian_day(const int day_const, const int month_const,
		  const int year_const)
{
    float a, b, jd;
    bool gregorian;

    int month = month_const;
    int day = day_const;
    float year = (float) year_const;

    gregorian = (year < 1583) ? 0 : 1;

    if ((month == 1) || (month == 2)) {
	year = year - 1;
	month = month + 12;
    }

    a = floor(year / 100);
    if (gregorian)
	b = 2.0 - a + floor(a / 4.0);
    else
	b = 0.0;

    jd = floor(365.25 * (float)(year + 4716.0))
	+ floor(30.6001 * (float)(month + 1))
	+ day + b - 1524.5;

    return jd;
}

// moon's position using fundamental arguments 
// (Van Flandern & Pulkkinen, 1979)
void moon(float jd)
{


  float d, f, g, h, m, n, s, u, v, w;

    h = 0.606434 + 0.03660110129 * jd;
    m = 0.374897 + 0.03629164709 * jd;
    f = 0.259091 + 0.0367481952 * jd;
    d = 0.827362 + 0.03386319198 * jd;
    n = 0.347343 - 0.00014709391 * jd;
    g = 0.993126 + 0.0027377785 * jd;

    h = h - floor(h);
    m = m - floor(m);
    f = f - floor(f);
    d = d - floor(d);
    n = n - floor(n);
    g = g - floor(g);

    h = h * 2 * M_PI;
    m = m * 2 * M_PI;
    f = f * 2 * M_PI;
    d = d * 2 * M_PI;
    n = n * 2 * M_PI;
    g = g * 2 * M_PI;

    v = 0.39558 * sin(f + n);
    v = v + 0.082 * sin(f);
    v = v + 0.03257 * sin(m - f - n);
    v = v + 0.01092 * sin(m + f + n);
    v = v + 0.00666 * sin(m - f);
    v = v - 0.00644 * sin(m + f - 2 * d + n);
    v = v - 0.00331 * sin(f - 2 * d + n);
    v = v - 0.00304 * sin(f - 2 * d);
    v = v - 0.0024 * sin(m - f - 2 * d - n);
    v = v + 0.00226 * sin(m + f);
    v = v - 0.00108 * sin(m + f - 2 * d);
    v = v - 0.00079 * sin(f - n);
    v = v + 0.00078 * sin(f + 2 * d + n);

    u = 1 - 0.10828 * cos(m);
    u = u - 0.0188 * cos(m - 2 * d);
    u = u - 0.01479 * cos(2 * d);
    u = u + 0.00181 * cos(2 * m - 2 * d);
    u = u - 0.00147 * cos(2 * m);
    u = u - 0.00105 * cos(2 * d - g);
    u = u - 0.00075 * cos(m - 2 * d + g);

    w = 0.10478 * sin(m);
    w = w - 0.04105 * sin(2 * f + 2 * n);
    w = w - 0.0213 * sin(m - 2 * d);
    w = w - 0.01779 * sin(2 * f + n);
    w = w + 0.01774 * sin(n);
    w = w + 0.00987 * sin(2 * d);
    w = w - 0.00338 * sin(m - 2 * f - 2 * n);
    w = w - 0.00309 * sin(g);
    w = w - 0.0019 * sin(2 * f);
    w = w - 0.00144 * sin(m + n);
    w = w - 0.00144 * sin(m - 2 * f - n);
    w = w - 0.00113 * sin(m + 2 * f + 2 * n);
    w = w - 0.00094 * sin(m - 2 * d + g);
    w = w - 0.00092 * sin(2 * m - 2 * d);

    s = w / sqrt(u - v * v);	// compute moon's right ascension ...  
    Sky[0] = h + atan(s / sqrt(1 - s * s));

    s = v / sqrt(u);		// declination ...
    Sky[1] = atan(s / sqrt(1 - s * s));

    Sky[2] = 60.40974 * sqrt(u);	// and parallax
}

// test an hour for an event
float test_moon(int k, float t0, float lat, float plx)
{

const static float K1 = 15 * M_PI * 1.0027379 / 180;
static float Rise_az = 0.0, Set_az = 0.0;

    float ha[3] = { 0.0, 0.0, 0.0 };
    float a, b, c, d, e, s, z;
    float hr, min, time;
    float az, hz, nz, dz;

    if (RAn[2] < RAn[0])
	RAn[2] = RAn[2] + 2 * M_PI;

    ha[0] = t0 - RAn[0] + (k * K1);
    ha[2] = t0 - RAn[2] + (k * K1) + K1;

    ha[1] = (ha[2] + ha[0]) / 2;	// hour angle at half hour
    Dec[1] = (Dec[2] + Dec[0]) / 2;	// declination at half hour

    s = sin(DR * lat);
    c = cos(DR * lat);

    // refraction + sun semidiameter at horizon + parallax correction
    z = cos(DR * (90.567 - 41.685 / plx));

    if (k <= 0)			// first call of function
	VHz[0] = s * sin(Dec[0]) + c * cos(Dec[0]) * cos(ha[0]) - z;

    VHz[2] = s * sin(Dec[2]) + c * cos(Dec[2]) * cos(ha[2]) - z;

    if (sgn(VHz[0]) == sgn(VHz[2]))
	return VHz[2];		// no event this hour

    VHz[1] = s * sin(Dec[1]) + c * cos(Dec[1]) * cos(ha[1]) - z;

    a = 2 * VHz[2] - 4 * VHz[1] + 2 * VHz[0];
    b = 4 * VHz[1] - 3 * VHz[0] - VHz[2];
    d = b * b - 4 * a * VHz[0];

    if (d < 0)
	return VHz[2];		// no event this hour

    d = sqrt(d);
    e = (-b + d) / (2 * a);

    if ((e > 1) || (e < 0))
	e = (-b - d) / (2 * a);

    time = ((float) k) + e + 1 / 120;	// time of an event + round up
    hr = floor(time);
    min = floor((time - hr) * 60);

    hz = ha[0] + e * (ha[2] - ha[0]);	// azimuth of the moon at the event
    nz = -cos(Dec[1]) * sin(hz);
    dz = c * sin(Dec[1]) - s * cos(Dec[1]) * cos(hz);
    az = atan2(nz, dz) / DR;
    if (az < 0)
	az = az + 360;

    if ((VHz[0] < 0) && (VHz[2] > 0)) {
	Rise_time[0] = (int) hr;
	Rise_time[1] = (int) min;
	Rise_az = az;
	Moonrise = 1;
    }

    if ((VHz[0] > 0) && (VHz[2] < 0)) {
	Set_time[0] = (int) hr;
	Set_time[1] = (int) min;
	Set_az = az;
	Moonset = 1;
    }

    return VHz[2];
}

// Local Sidereal Time for zone
float lst(const float lon, const float jd, const float z)
{
    float s =
	24110.5 + 8640184.812999999 * jd / 36525 + 86636.6 * z +
	86400 * lon;
    s = s / 86400;
    s = s - floor(s);
    return s * 360 * DR;
}

// 3-point interpolation
float interpolate(const float f0, const float f1, const float f2,
		   const float p)
{
    float a = f1 - f0;
    float b = f2 - f1 - a;
    float f = f0 + p * (2 * a + b * (2 * p - 1));

    return f;
}


// calculate moonrise and moonset times
void riseset(const float lat, const float lon, const int day,
	     const int month, const int year, const int TimezoneOffset)
{
    int i, j, k;
    float ph;
    // guido: julian day has been converted to int from float
//    float jd = (julian_day(day, month, year)) - 2451545;	// Julian day relative to Jan 1.5, 2000
    float jd = (julian_day(day, month, year)) - 2451545;	// Julian day relative to Jan 1.5, 2000
    float mp[3][3];
    float lon_local = lon;

    if ((sgn(-TimezoneOffset) == sgn(lon)) && (TimezoneOffset != 0))
	Serial.println("WARNING: time zone and longitude are incompatible!");

    for (i = 0; i < 3; i++) {
	for (j = 0; j < 3; j++)
	    mp[i][j] = 0.0;
    }

    lon_local = lon / 360;
    float tz = -((float)TimezoneOffset) / 24;
    float t0 = lst(lon_local, jd, tz);	// local sidereal time

    jd = jd + tz;		// get moon position at start of day

    for (k = 0; k < 3; k++) {
	moon(jd);
	mp[k][0] = Sky[0];
	mp[k][1] = Sky[1];
	mp[k][2] = Sky[2];
	jd = jd + 0.5;
    }

    if (mp[1][0] <= mp[0][0])
	mp[1][0] = mp[1][0] + 2 * M_PI;

    if (mp[2][0] <= mp[1][0])
	mp[2][0] = mp[2][0] + 2 * M_PI;

    RAn[0] = mp[0][0];
    Dec[0] = mp[0][1];

    Moonrise = 0;		// initialize
    Moonset = 0;

    for (k = 0; k < 24; k++)	// check each hour of this day
    {
	ph = ((float) (k + 1)) / 24;

	RAn[2] = interpolate(mp[0][0], mp[1][0], mp[2][0], ph);
	Dec[2] = interpolate(mp[0][1], mp[1][1], mp[2][1], ph);

	VHz[2] = test_moon(k, t0, lat, mp[1][2]);

	RAn[0] = RAn[2];	// advance to next hour
	Dec[0] = Dec[2];
	VHz[0] = VHz[2];
    }

}

float get_phase(const DateTime now, const int TZOffset) {
  //float phase = (julian_day(now.day(),now.month(), now.year()) - 2451550.1);
  //float phase = (julian_day(now.day(),now.month(), now.year()) - 2440594.359028);  //CE 1970 January 07 20:37:00.2 UT  Wednesday - new moon
  float phase;
  
  phase = julian_day(now.day(),now.month(), now.year());
//  phase -= 2455211.8; //CE  2010 January 15 07:12:00.0 UT - new moon
  phase -= 2440594.359028;
  phase += ((now.hour() * 60) + now.minute()) / (float)1440; // Adjust for current time JD returns Midnight. 1140 = minutes in day.
  phase += -((float)TZOffset) / 24.0; // Now adjust local to UTC.
  phase = fmod(phase,29.530588853);
  return (phase);
}


// check for no moonrise and/or no moonset
int moon_up(const DateTime now) {
  
  int riseMin=(Rise_time[0]*60)+Rise_time[1];
  int setMin=(Set_time[0]*60)+Set_time[1];
  int nowMin=(now.hour()*60)+now.minute();

  if ((!Moonrise) && (!Moonset)) { // neither moonrise nor moonset
    if (VHz[2] < 0)
      return(0); // down all day
    else
      return(1); // up all day
    }
  
  if (Moonrise && Moonset) {
    if ((setMin > riseMin) && (riseMin < nowMin) && (nowMin < setMin)) {
      return(1); // up
    } 

    if ((setMin < riseMin) && ((nowMin < setMin) || (nowMin > riseMin))) {
      return(1); // up
    } 
  }
  
  if (Moonrise && (!Moonset)) { // Moonrise only
    if (nowMin > riseMin) {
      return(1);
    }
  }
  
  if (Moonset && (!Moonrise)) { // Moonset only
    if (nowMin < setMin) {
      return(1);
    }
  }

  return(0); // if in doubt turn blow it out

}


void print_time(const DateTime now) {
  Serial.print(now.year(), DEC);
  Serial.print('/');
  if (now.year() < 10)
    Serial.print("0");
  Serial.print(now.month(), DEC);
  Serial.print('/');
  if (now.day() < 10)
    Serial.print("0");
  Serial.print(now.day(), DEC);
  Serial.print(' ');
  if (now.hour() < 10)
    Serial.print("0");
  Serial.print(now.hour(), DEC);
  Serial.print(':');
  if (now.minute() < 10)
    Serial.print("0");
  Serial.print(now.minute(), DEC);
  Serial.print(':');
  if (now.second() < 10)
    Serial.print("0");
  Serial.print(now.second(), DEC);
}

void print_moonrise (void) {
  if (Moonrise) {
    if (Rise_time[0] < 10) {
      Serial.print("0");
    }
    Serial.print(Rise_time[0]);
    Serial.print(":");
    if (Rise_time[1] < 10) {
      Serial.print("0");
    }
    Serial.print(Rise_time[1]);
  } else {
    Serial.print("  none");
  }
}

void print_moonset (void) {
  if (Moonset) {
    if (Set_time[0] < 10) {
      Serial.print("0");
    }
    Serial.print(Set_time[0]);
    Serial.print(":");
    if (Set_time[1] < 10) {
      Serial.print("0");
    }
    Serial.print(Set_time[1]);
  } else {
    Serial.print("  none");
  }
}

void print_phase(const DateTime now) {
  Serial.print(get_phase(now, TZOFFSET),4);
}


void show_phase(const float phase) {
  int phase_int = int(phase*(12/29.53));
  int phase_pwm = int(((phase*(12/29.53))-phase_int) * float(LEDMAX));

  switch(phase_int) {
    case 0:
      analogWrite(LED1, phase_pwm);
      analogWrite(LED2, 0);
      analogWrite(LED3, 0);
      analogWrite(LED4, 0);
      analogWrite(LED5, 0);
      analogWrite(LED6, 0);
      break;
    case 1:
      analogWrite(LED1, LEDMAX);
      analogWrite(LED2, phase_pwm);
      analogWrite(LED3, 0);
      analogWrite(LED4, 0);
      analogWrite(LED5, 0);
      analogWrite(LED6, 0);
      break;
    case 2:
      analogWrite(LED1, LEDMAX);
      analogWrite(LED2, LEDMAX);
      analogWrite(LED3, phase_pwm);
      analogWrite(LED4, 0);
      analogWrite(LED5, 0);
      analogWrite(LED6, 0);
      break;
    case 3:
      analogWrite(LED1, LEDMAX);
      analogWrite(LED2, LEDMAX);
      analogWrite(LED3, LEDMAX);
      analogWrite(LED4, phase_pwm);
      analogWrite(LED5, 0);
      analogWrite(LED6, 0);
      break;
    case 4:

      analogWrite(LED1, LEDMAX);
      analogWrite(LED2, LEDMAX);
      analogWrite(LED3, LEDMAX);
      analogWrite(LED4, LEDMAX);
      analogWrite(LED5, phase_pwm);
      analogWrite(LED6, 0);
      break;
    case 5:
      analogWrite(LED1, LEDMAX);
      analogWrite(LED2, LEDMAX);
      analogWrite(LED3, LEDMAX);
      analogWrite(LED4, LEDMAX);
      analogWrite(LED5, LEDMAX);
      analogWrite(LED6, phase_pwm);
      break;
    case 6:
      analogWrite(LED1, LEDMAX-phase_pwm);
      analogWrite(LED2, LEDMAX);
      analogWrite(LED3, LEDMAX);
      analogWrite(LED4, LEDMAX);
      analogWrite(LED5, LEDMAX);
      analogWrite(LED6, LEDMAX);
      break;
    case 7:
      analogWrite(LED1, 0);
      analogWrite(LED2, LEDMAX-phase_pwm);
      analogWrite(LED3, LEDMAX);
      analogWrite(LED4, LEDMAX);
      analogWrite(LED5, LEDMAX);
      analogWrite(LED6, LEDMAX);
      break;
    case 8:
      analogWrite(LED1, 0);
      analogWrite(LED2, 0);
      analogWrite(LED3, LEDMAX-phase_pwm);
      analogWrite(LED4, LEDMAX);
      analogWrite(LED5, LEDMAX);
      analogWrite(LED6, LEDMAX);
      break;
    case 9:
      analogWrite(LED1, 0);
      analogWrite(LED2, 0);
      analogWrite(LED3, 0);
      analogWrite(LED4, LEDMAX-phase_pwm);
      analogWrite(LED5, LEDMAX);
      analogWrite(LED6, LEDMAX);
      break;
    case 10:
      analogWrite(LED1, 0);
      analogWrite(LED2, 0);
      analogWrite(LED3, 0);
      analogWrite(LED4, 0);
      analogWrite(LED5, LEDMAX-phase_pwm);
      analogWrite(LED6, LEDMAX);
      break;
    case 11:
      analogWrite(LED1, 0);
      analogWrite(LED2, 0);
      analogWrite(LED3, 0);
      analogWrite(LED4, 0);
      analogWrite(LED5, 0);
      analogWrite(LED6, LEDMAX-phase_pwm);
      break;    
  }
  
}
void show_demo(void) {

  int x,y;

  for(y=1000; y > 400; y-=200) {
    for (x=0;x<=2953;x++) {
      show_phase(x/100.0); 
      delayMicroseconds(y);
    }
  }

  for(y=200; y > 0; y-=25) {
    for (x=0;x<=2953;x++) {
      show_phase(x/100.0); 
      delayMicroseconds(y);
    }
  }

  for(y=6000;y>1000;y-=1000) {
    for(x=0; x<=LEDMAX; x++) {
      analogWrite(LED1, x);  
      analogWrite(LED2, x);  
      analogWrite(LED3, x);  
      analogWrite(LED4, x);  
      analogWrite(LED5, x);  
      analogWrite(LED6, x);  
      delayMicroseconds(y);
    }
    delay(5);
    for(x=LEDMAX; x>=0; x--) {
      analogWrite(LED1, x);  
      analogWrite(LED2, x);  
      analogWrite(LED3, x);  
      analogWrite(LED4, x);  
      analogWrite(LED5, x);  
      analogWrite(LED6, x);  
      delayMicroseconds(y);
    }
    delay(5);
  }

  for(y=1000;y<6000;y+=1000) {
    for(x=0; x<=LEDMAX; x++) {
      analogWrite(LED1, x);  
      analogWrite(LED2, x);  
      analogWrite(LED3, x);  
      analogWrite(LED4, x);  
      analogWrite(LED5, x);  
      analogWrite(LED6, x);  
      delayMicroseconds(y);
    }
    delay(5);
    for(x=255; x>=0; x--) {
      analogWrite(LED1, x);  
      analogWrite(LED2, x);  
      analogWrite(LED3, x);  
      analogWrite(LED4, x);  
      analogWrite(LED5, x);  
      analogWrite(LED6, x);  
      delayMicroseconds(y);
    }
    delay(5);
  }
  delay(1000);
}

void show_error(void) {
  int x;
  digitalWrite(LED1, HIGH);  
  digitalWrite(LED2, HIGH);  
  digitalWrite(LED3, HIGH);  
  digitalWrite(LED4, HIGH);  
  digitalWrite(LED5, HIGH);  
  digitalWrite(LED6, HIGH);  
  delay(500);
  digitalWrite(LED1, LOW);  
  digitalWrite(LED2, LOW);  
  digitalWrite(LED3, LOW);  
  digitalWrite(LED4, LOW);  
  digitalWrite(LED5, LOW);  
  digitalWrite(LED6, LOW);
  delay(500);  
}

void setup () {
  Serial.begin(57600);
  Wire.begin();
  RTC.begin();
    
  pinMode(LED1, OUTPUT); 
  pinMode(LED2, OUTPUT); 
  pinMode(LED3, OUTPUT); 
  pinMode(LED4, OUTPUT); 
  pinMode(LED5, OUTPUT); 
  pinMode(LED6, OUTPUT); 
}

void loop () {
  if (RTC.isrunning()) {
    DateTime now = RTC.now();

    if (start) {
      show_demo();
      start=0;
    }
    if (then != now.second()) {
  
      riseset(LATITUDE, LONGITUDE, now.day(), now.month(), now.year(), TZOFFSET);
  
      if(moon_up(now)) {
        show_phase(get_phase(now, TZOFFSET));
      } else {
        show_phase(0);
      }
  
      Serial.print(" Date/Time: ");
      print_time(now);
      Serial.print(" Moonrise: ");
      print_moonrise();
      Serial.print(" Moonset: ");
      print_moonset();
      Serial.print(" Age: ");
      print_phase(now);
      Serial.print(" State: ");
      Serial.print(moon_up(now)?"Up":"Down");
      Serial.print("\n");
      then=now.second();
    }
  } else {
    show_error();
  }
}

