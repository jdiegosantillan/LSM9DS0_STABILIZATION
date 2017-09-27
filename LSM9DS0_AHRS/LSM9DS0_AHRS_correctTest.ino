/*****************************************************************
LSM9DS0_AHRS.ino
SFE_LSM9DS0 Library AHRS Data Fusion Example Code
Jim Lindblom @ SparkFun Electronics
Original Creation Date: February 18, 2014
https://github.com/sparkfun/LSM9DS0_Breakout

Modified by Kris Winer, April 4, 2014

The LSM9DS0 is a versatile 9DOF sensor. It has a built-in
accelerometer, gyroscope, and magnetometer. Very cool! Plus it
functions over either SPI or I2C.

This Arduino sketch utilizes Jim Lindblom's SFE_LSM9DS0 library to generate the basic sensor data
for use in two sensor fusion algorithms becoming increasingly popular with DIY quadcopter and robotics engineers. 

Like the original LSM9SD0_simple.ino sketch, it'll demo the following:
* How to create a LSM9DS0 object, using a constructor (global
  variables section).
* How to use the begin() function of the LSM9DS0 class.
* How to read the gyroscope, accelerometer, and magnetometer
  using the readGryo(), readAccel(), readMag() functions and the
  gx, gy, gz, ax, ay, az, mx, my, and mz variables.
* How to calculate actual acceleration, rotation speed, magnetic
  field strength using the calcAccel(), calcGyro() and calcMag()
  functions.
  
In addition, the sketch will demo:
* How to check for data updates using interrupts
* How to display output at a rate different from the sensor data update and fusion filter update rates
* How to specify the accelerometer anti-aliasing (low-pass) filter rate
* How to use the data from the LSM9DS0 to fuse the sensor data into a quaternion representation of the sensor frame
  orientation relative to a fixed Earth frame providing absolute orientation information for subsequent use.
* An example of how to use the quaternion data to generate standard aircraft orientation data in the form of
  Tait-Bryan angles representing the sensor yaw, pitch, and roll angles suitable for any vehicle stablization control application.

Hardware setup: This library supports communicating with the
LSM9DS0 over either I2C or SPI. If you're using I2C, these are
the only connections that need to be made:
	LSM9DS0 --------- Arduino
	 SCL ---------- SCL (A5 on older 'Duinos')
	 SDA ---------- SDA (A4 on older 'Duinos')
	 VDD ------------- 3.3V
	 GND ------------- GND
         DRDYG-------------49   (gyro data ready interrupt, can be any digital pin)
         INTX1-------------53   (accelerometer data ready interrupt, can be any digital pin)
         INTX2-------------51   (magnetometer data ready interrupt, can be any digital pin)
(CSG, CSXM, SDOG, and SDOXM should all be pulled high jumpers on 
  the breakout board will do this for you.) In other words put these pins to 3.3 V 
  
 Note: The LSM9DS0 in the I2C mode uses the Arduino Wire library. 
 Because the sensor is not 5V tolerant, we are using a 3.3 V 8 MHz Pro Mini or a 3.3 V Teensy 3.1.
 We have disabled the internal pull-ups used by the Wire library in the Wire.h/twi.c utility file.
 We are also using the 400 kHz fast I2C mode by setting the TWI_FREQ  to 400000L /twi.h utility file.
  

This code is beerware. If you see me (or any other SparkFun 
employee) at the local, and you've found our code helpful, please 
buy us a round!

Distributed as-is; no warranty is given.
*****************************************************************/

// The SFE_LSM9DS0 requires both the SPI and Wire libraries.
// Unfortunately, you'll need to include both in the Arduino
// sketch, before including the SFE_LSM9DS0 library.
#include <SPI.h> // Included for SFE_LSM9DS0 library
#include <Wire.h>
#include "SFE_LSM9DS0.h"


///////////////////////
// Example I2C Setup //
///////////////////////
// Comment out this section if you're using SPI
// SDO_XM and SDO_G are both grounded, so our addresses are:
#define LSM9DS0_XM  0x1D // Would be 0x1E if SDO_XM is LOW
#define LSM9DS0_G   0x6B // Would be 0x6A if SDO_G is LOW
// Create an instance of the LSM9DS0 library called `dof` the
// parameters for this constructor are:
// [SPI or I2C Mode declaration],[gyro I2C address],[xm I2C add.]
LSM9DS0 dof(MODE_I2C, LSM9DS0_G, LSM9DS0_XM);


///////////////////////////////
// Interrupt Pin Definitions //
///////////////////////////////
const byte INT1XM = 53; // INT1XM tells us when accel data is ready
const byte INT2XM = 51; // INT2XM tells us when mag data is ready
const byte DRDYG  = 49; // DRDYG  tells us when gyro data is ready

// global constants for 9 DoF fusion and AHRS (Attitude and Heading Reference System)
#define GyroMeasError PI * (40.0f / 180.0f)       // gyroscope measurement error in rads/s (shown as 3 deg/s)
#define GyroMeasDrift PI * (0.0f / 180.0f)      // gyroscope measurement drift in rad/s/s (shown as 0.0 deg/s/s)
// There is a tradeoff in the beta parameter between accuracy and response speed.
// In the original Madgwick study, beta of 0.041 (corresponding to GyroMeasError of 2.7 degrees/s) was found to give optimal accuracy.
// However, with this value, the LSM9SD0 response time is about 10 seconds to a stable initial quaternion.
// Subsequent changes also require a longish lag time to a stable output, not fast enough for a quadcopter or robot car!
// By increasing beta (GyroMeasError) by about a factor of fifteen, the response time constant is reduced to ~2 sec
// I haven't noticed any reduction in solution accuracy. This is essentially the I coefficient in a PID control sense; 
// the bigger the feedback coefficient, the faster the solution converges, usually at the expense of accuracy. 
// In any case, this is the free parameter in the Madgwick filtering and fusion scheme.
#define beta sqrt(3.0f / 4.0f) * GyroMeasError   // compute beta
#define zeta sqrt(3.0f / 4.0f) * GyroMeasDrift   // compute zeta, the other free parameter in the Madgwick scheme usually set to a small or zero value
#define Kp 2.0f * 5.0f // these are the free parameters in the Mahony filter and fusion scheme, Kp for proportional feedback, Ki for integral
#define Ki 0.0f

uint32_t count = 0;  // used to control display output rate
uint32_t delt_t = 0; // used to control display output rate
float pitch, yaw, roll, heading;
float deltat = 0.0f;        // integration interval for both filter schemes
uint32_t lastUpdate = 0;    // used to calculate integration interval
uint32_t Now = 0;           // used to calculate integration interval

float abias[3] = {0, 0, 0}, gbias[3] = {0, 0, 0};
float ax, ay, az, gx, gy, gz, mx, my, mz; // variables to hold latest sensor data values 
float q[4] = {1.0f, 0.0f, 0.0f, 0.0f};    // vector to hold quaternion
float eInt[3] = {0.0f, 0.0f, 0.0f};       // vector to hold integral error for Mahony method
float temperature;
int i; 

void setup()
{
  Serial.begin(38400); // Start serial at 38400 bps
 
  // Set up interrupt pins as inputs:
  pinMode(INT1XM, INPUT);
  pinMode(INT2XM, INPUT);
  pinMode(DRDYG,  INPUT);

  pinMode(2, OUTPUT); //1 on motor  
  pinMode(3, OUTPUT); //2 on motor
  pinMode(4, OUTPUT); //3 on motor
  pinMode(5, OUTPUT); //4 on motor
  pinMode(6, OUTPUT); //5 on motor
  pinMode(7, OUTPUT); //6 on motor
  pinMode(8, OUTPUT); //7 on motor
  pinMode(9, OUTPUT); //8 on motor

  
  // begin() returns a 16-bit value which includes both the gyro 
  // and accelerometers WHO_AM_I response. You can check this to
  // make sure communication was successful.
  uint32_t status = dof.begin();
  
 // Set data output ranges; choose lowest ranges for maximum resolution
 // Accelerometer scale can be: A_SCALE_2G, A_SCALE_4G, A_SCALE_6G, A_SCALE_8G, or A_SCALE_16G   
    dof.setAccelScale(dof.A_SCALE_2G);
 // Gyro scale can be:  G_SCALE__245, G_SCALE__500, or G_SCALE__2000DPS
    dof.setGyroScale(dof.G_SCALE_245DPS);
 // Magnetometer scale can be: M_SCALE_2GS, M_SCALE_4GS, M_SCALE_8GS, M_SCALE_12GS   
    dof.setMagScale(dof.M_SCALE_2GS);
  
 // Set output data rates  
 // Accelerometer output data rate (ODR) can be: A_ODR_3125 (3.225 Hz), A_ODR_625 (6.25 Hz), A_ODR_125 (12.5 Hz), A_ODR_25, A_ODR_50, 
 //                                              A_ODR_100,  A_ODR_200, A_ODR_400, A_ODR_800, A_ODR_1600 (1600 Hz)
    dof.setAccelODR(dof.A_ODR_200); // Set accelerometer update rate at 100 Hz
 // Accelerometer anti-aliasing filter rate can be 50, 194, 362, or 763 Hz
 // Anti-aliasing acts like a low-pass filter allowing oversampling of accelerometer and rejection of high-frequency spurious noise.
 // Strategy here is to effectively oversample accelerometer at 100 Hz and use a 50 Hz anti-aliasing (low-pass) filter frequency
 // to get a smooth ~150 Hz filter update rate
    dof.setAccelABW(dof.A_ABW_50); // Choose lowest filter setting for low noise
 
 // Gyro output data rates can be: 95 Hz (bandwidth 12.5 or 25 Hz), 190 Hz (bandwidth 12.5, 25, 50, or 70 Hz)
 //                                 380 Hz (bandwidth 20, 25, 50, 100 Hz), or 760 Hz (bandwidth 30, 35, 50, 100 Hz)
    dof.setGyroODR(dof.G_ODR_190_BW_125);  // Set gyro update rate to 190 Hz with the smallest bandwidth for low noise

 // Magnetometer output data rate can be: 3.125 (ODR_3125), 6.25 (ODR_625), 12.5 (ODR_125), 25, 50, or 100 Hz
    dof.setMagODR(dof.M_ODR_125); // Set magnetometer to update every 80 ms
    
 // Use the FIFO mode to average accelerometer and gyro readings to calculate the biases, which can then be removed from
 // all subsequent measurements.
    dof.calLSM9DS0(gbias, abias);
}

void loop()
{
  if(digitalRead(DRDYG)) {  // When new gyro data is ready
  dof.readGyro();           // Read raw gyro data
    gx = dof.calcGyro(dof.gx) - gbias[0];   // Convert to degrees per seconds, remove gyro biases
    gy = dof.calcGyro(dof.gy) - gbias[1];
    gz = dof.calcGyro(dof.gz) - gbias[2];
  }
  
  if(digitalRead(INT1XM)) {  // When new accelerometer data is ready
    dof.readAccel();         // Read raw accelerometer data
    ax = dof.calcAccel(dof.ax) - abias[0];   // Convert to g's, remove accelerometer biases
    ay = dof.calcAccel(dof.ay) - abias[1];
    az = dof.calcAccel(dof.az) - abias[2];
  }
  
  if(digitalRead(INT2XM)) {  // When new magnetometer data is ready
    dof.readMag();           // Read raw magnetometer data
    mx = dof.calcMag(dof.mx);     // Convert to Gauss and correct for calibration
    my = dof.calcMag(dof.my);
    mz = dof.calcMag(dof.mz);
    
    dof.readTemp();
    temperature = 21.0 + (float) dof.temperature/8.; // slope is 8 LSB per degree C, just guessing at the intercept
  }

  Now = micros();
  deltat = ((Now - lastUpdate)/1000000.0f); // set integration time by time elapsed since last filter update
  lastUpdate = Now;
  // Sensors x- and y-axes are aligned but magnetometer z-axis (+ down) is opposite to z-axis (+ up) of accelerometer and gyro!
  // This is ok by aircraft orientation standards!  
  // Pass gyro rate as rad/s

 MahonyQuaternionUpdate(ax, ay, az, gx*PI/180.0f, gy*PI/180.0f, gz*PI/180.0f, mx, my, mz);

    // Serial print and/or display at 0.5 s rate independent of data rates
    delt_t = millis() - count;
    if (delt_t > 10) { // update LCD once per half-second independent of read rate
 

  // Define output variables from updated quaternion---these are Tait-Bryan angles, commonly used in aircraft orientation.
  // In this coordinate system, the positive z-axis is down toward Earth. 
  // Yaw is the angle between Sensor x-axis and Earth magnetic North (or true North if corrected for local declination), 
  // looking down on the sensor positive yaw is counterclockwise.
  // Pitch is angle between sensor x-axis and Earth ground plane, toward the Earth is positive, up toward the sky is negative.
  // Roll is angle between sensor y-axis and Earth ground plane, y-axis up is positive roll.
  // These arise from the definition of the homogeneous rotation matrix constructed from quaternions.
  // Tait-Bryan angles as well as Euler angles are non-commutative; that is, to get the correct orientation the rotations must be
  // applied in the correct order which for this configuration is yaw, pitch, and then roll.
  // For more see http://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles which has additional links.
    yaw   = atan2(2.0f * (q[1] * q[2] + q[0] * q[3]), q[0] * q[0] + q[1] * q[1] - q[2] * q[2] - q[3] * q[3]);   
    pitch = -asin(2.0f * (q[1] * q[3] - q[0] * q[2]));
    roll  = atan2(2.0f * (q[0] * q[1] + q[2] * q[3]), q[0] * q[0] - q[1] * q[1] - q[2] * q[2] + q[3] * q[3]);
    pitch *= 180.0f / PI;
    yaw   *= 180.0f / PI; 
    yaw   -= 13.8; // Declination at Danville, California is 13 degrees 48 minutes and 47 seconds on 2014-04-04
    roll  *= 180.0f / PI;

    pitch+= - 7.3;

    Serial.print(roll); Serial.print('\t'); Serial.print(pitch); Serial.print('\t'); Serial.println(yaw);
  
    // With ODR settings of 400 Hz, 380 Hz, and 25 Hz for the accelerometer, gyro, and magnetometer, respectively,
    // the filter is updating at a ~125 Hz rate using the Madgwick scheme and ~165 Hz using the Mahony scheme 
    // even though the display refreshes at only 2 Hz.
    // The filter update rate can be increased by reducing the rate of data reading. The optimal implementation is
    // one which balances the competing rates so they are about the same; that is, the filter updates the sensor orientation
    // at about the same rate the data is changing. Of course, other implementations are possible. One might consider
    // updating the filter at twice the average new data rate to allow for finite filter convergence times.
    // The filter update rate is determined mostly by the mathematical steps in the respective algorithms, 
    // the processor speed (8 MHz for the 3.3V Pro Mini), and the sensor ODRs, especially the magnetometer ODR:
    // smaller ODRs for the magnetometer produce the above rates, maximum magnetometer ODR of 100 Hz produces
    // filter update rates of ~110 and ~135 Hz for the Madgwick and Mahony schemes, respectively. 
    // This is presumably because the magnetometer read takes longer than the gyro or accelerometer reads.
    // With low ODR settings of 100 Hz, 95 Hz, and 6.25 Hz for the accelerometer, gyro, and magnetometer, respectively,
    // the filter is updating at a ~150 Hz rate using the Madgwick scheme and ~200 Hz using the Mahony scheme.
    // These filter update rates should be fast enough to maintain accurate platform orientation for 
    // stabilization control of a fast-moving robot or quadcopter. Compare to the update rate of 200 Hz
    // produced by the on-board Digital Motion Processor of Invensense's MPU6050 6 DoF and MPU9150 9DoF sensors.
    // The 3.3 V 8 MHz Pro Mini is doing pretty well!
    stayLeveling();
    count = millis();
    }
}  
 // Similar to Madgwick scheme but uses proportional and integral filtering on the error between estimated reference vectors and
 // measured ones. 
 void stayLeveling(){
  for (i = 0; (2 * i) < 90; i++){ //loop will start from 0 degrees -> 90 degrees 
    //left(positive)
    if((roll > 2*i) && (roll < (2*i + 2))){
      //Boost the right motors
     analogWrite(4, i*5);
     analogWrite(5, i*5);
     analogWrite(6, i*2);
     analogWrite(7, i*2);
    }
    //right(negative)
    if((roll < -1 *(2*i)) && (roll > -1 *(2*i + 2))){
      //Boost the left motors
     analogWrite(2, i*5);
     analogWrite(3, i*5);
     analogWrite(8, i*2);
     analogWrite(9, i*2);  
    }
    //backward(positive)
    if((pitch > 2*i) && (pitch < (2*i + 2))){
      //Boost the back motors
     analogWrite(3, i*5);
     analogWrite(5, i*5);
     analogWrite(6, i*2);
     analogWrite(8, i*2);    
    }
    //forward
    if((pitch < -1*( 2*i)) && (pitch > -1 *(2*i + 2))){
      //Boost the front motors
     analogWrite(2, i*5);
     analogWrite(4, i*5);
     analogWrite(7, i*2);
     analogWrite(9, i*2);
    }
  }
 }
            void MahonyQuaternionUpdate(float ax, float ay, float az, float gx, float gy, float gz, float mx, float my, float mz)
        {
            float q1 = q[0], q2 = q[1], q3 = q[2], q4 = q[3];   // short name local variable for readability
            float norm;
            float hx, hy, bx, bz;
            float vx, vy, vz, wx, wy, wz;
            float ex, ey, ez;
            float pa, pb, pc;

            // Auxiliary variables to avoid repeated arithmetic
            float q1q1 = q1 * q1;
            float q1q2 = q1 * q2;
            float q1q3 = q1 * q3;
            float q1q4 = q1 * q4;
            float q2q2 = q2 * q2;
            float q2q3 = q2 * q3;
            float q2q4 = q2 * q4;
            float q3q3 = q3 * q3;
            float q3q4 = q3 * q4;
            float q4q4 = q4 * q4;   

            // Normalise accelerometer measurement
            norm = sqrt(ax * ax + ay * ay + az * az);
            if (norm == 0.0f) return; // handle NaN
            norm = 1.0f / norm;        // use reciprocal for division
            ax *= norm;
            ay *= norm;
            az *= norm;

            // Normalise magnetometer measurement
            norm = sqrt(mx * mx + my * my + mz * mz);
            if (norm == 0.0f) return; // handle NaN
            norm = 1.0f / norm;        // use reciprocal for division
            mx *= norm;
            my *= norm;
            mz *= norm;

            // Reference direction of Earth's magnetic field
            hx = 2.0f * mx * (0.5f - q3q3 - q4q4) + 2.0f * my * (q2q3 - q1q4) + 2.0f * mz * (q2q4 + q1q3);
            hy = 2.0f * mx * (q2q3 + q1q4) + 2.0f * my * (0.5f - q2q2 - q4q4) + 2.0f * mz * (q3q4 - q1q2);
            bx = sqrt((hx * hx) + (hy * hy));
            bz = 2.0f * mx * (q2q4 - q1q3) + 2.0f * my * (q3q4 + q1q2) + 2.0f * mz * (0.5f - q2q2 - q3q3);

            // Estimated direction of gravity and magnetic field
            vx = 2.0f * (q2q4 - q1q3);
            vy = 2.0f * (q1q2 + q3q4);
            vz = q1q1 - q2q2 - q3q3 + q4q4;
            wx = 2.0f * bx * (0.5f - q3q3 - q4q4) + 2.0f * bz * (q2q4 - q1q3);
            wy = 2.0f * bx * (q2q3 - q1q4) + 2.0f * bz * (q1q2 + q3q4);
            wz = 2.0f * bx * (q1q3 + q2q4) + 2.0f * bz * (0.5f - q2q2 - q3q3);  

            // Error is cross product between estimated direction and measured direction of gravity
            ex = (ay * vz - az * vy) + (my * wz - mz * wy);
            ey = (az * vx - ax * vz) + (mz * wx - mx * wz);
            ez = (ax * vy - ay * vx) + (mx * wy - my * wx);
            if (Ki > 0.0f)
            {
                eInt[0] += ex;      // accumulate integral error
                eInt[1] += ey;
                eInt[2] += ez;
            }
            else
            {
                eInt[0] = 0.0f;     // prevent integral wind up
                eInt[1] = 0.0f;
                eInt[2] = 0.0f;
            }

            // Apply feedback terms
            gx = gx + Kp * ex + Ki * eInt[0];
            gy = gy + Kp * ey + Ki * eInt[1];
            gz = gz + Kp * ez + Ki * eInt[2];

            // Integrate rate of change of quaternion
            pa = q2;
            pb = q3;
            pc = q4;
            q1 = q1 + (-q2 * gx - q3 * gy - q4 * gz) * (0.5f * deltat);
            q2 = pa + (q1 * gx + pb * gz - pc * gy) * (0.5f * deltat);
            q3 = pb + (q1 * gy - pa * gz + pc * gx) * (0.5f * deltat);
            q4 = pc + (q1 * gz + pa * gy - pb * gx) * (0.5f * deltat);

            // Normalise quaternion
            norm = sqrt(q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4);
            norm = 1.0f / norm;
            q[0] = q1 * norm;
            q[1] = q2 * norm;
            q[2] = q3 * norm;
            q[3] = q4 * norm;
 
        }

