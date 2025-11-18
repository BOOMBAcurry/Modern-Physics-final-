/*
 * LED Planck Constant Measurement System
 * 
 * Hardware Setup:
 * - MCP4725 DAC: Controls LED voltage (I2C address 0x60 or 0x61)
 * - AS7341 Spectrometer: Measures LED wavelength (I2C address 0x39)
 * - LED Circuit: DAC -> Resistor (1200Ω) -> A0 pin -> LED Anode -> LED Cathode -> GND
 * 
 * The A0 pin measures voltage at the LED anode to calculate current through the LED
 * 
 * I2C Connections (both devices share these pins):
 * - SDA: Arduino A4 (Uno/Nano) or pin 20 (Mega)
 * - SCL: Arduino A5 (Uno/Nano) or pin 21 (Mega)
 */

#include <Wire.h>
#include <Adafruit_MCP4725.h>
#include <Adafruit_AS7341.h>

// Create objects for DAC and spectrometer
Adafruit_MCP4725 dac;
Adafruit_AS7341 as;

// ===== PIN DEFINITIONS =====
const int measure_led_v = A0;      // Analog pin to measure LED voltage

// ===== CIRCUIT CONSTANTS =====
const float Vcc_nominal = 5.0;     // Arduino supply voltage (change to 3.3 if using 3.3V)
const float ADC_COUNTS = 1023.0;   // 10-bit ADC (0-1023)
const int DAC_COUNTS = 4095;       // 12-bit DAC (0-4095)
const float R_series = 1200.0;     // Series resistor in ohms (CHANGE THIS TO YOUR RESISTOR VALUE)

// ===== SWEEP PARAMETERS =====
const int sweep_step = 8;          // DAC increment per step (larger = faster but less detailed)
const int settle_ms = 5;           // Wait time after changing DAC for voltage to stabilize

// ===== WAVELENGTH CHANNEL CENTERS (nm) =====
// AS7341 has 8 channels (F1-F8) centered at these wavelengths
const float channel_centers_nm[8] = {
  415.0,  // F1 - Violet
  445.0,  // F2 - Blue
  480.0,  // F3 - Cyan
  515.0,  // F4 - Green
  555.0,  // F5 - Yellow-Green
  590.0,  // F6 - Yellow
  630.0,  // F7 - Orange-Red
  680.0   // F8 - Red
};

// ===== GLOBAL VARIABLES =====
bool correct_Vcc = false;          // Enable/disable Vcc correction
String serial_command;              // Stores incoming serial commands

// ===== HELPER FUNCTIONS =====

/**
 * Read analog voltage from a pin
 * @param pin - Analog pin number (e.g., A0)
 * @param vcc - Reference voltage for ADC
 * @return Voltage in volts
 */
float readAnalogVoltage(int pin, float vcc) {
  int raw = analogRead(pin);
  return (raw / ADC_COUNTS) * vcc;
}

/**
 * Read actual Vcc voltage using internal bandgap reference
 * This compensates for USB voltage variations
 * @return Vcc in millivolts
 */
long readVcc() {
  long result;
  // Set ADC to read internal 1.1V reference against Vcc
#if defined(__AVR_ATmega32U4__) || defined(__AVR_ATmega1280__) || defined(__AVR_ATmega2560__)
  ADMUX = _BV(REFS0) | _BV(MUX4) | _BV(MUX3) | _BV(MUX2) | _BV(MUX1);
#elif defined (__AVR_ATtiny24__) || defined(__AVR_ATtiny44__) || defined(__AVR_ATtiny84__)
  ADMUX = _BV(MUX5) | _BV(MUX0);
#elif defined (__AVR_ATtiny25__) || defined(__AVR_ATtiny45__) || defined(__AVR_ATtiny85__)
  ADMUX = _BV(MUX3) | _BV(MUX2);
#else
  ADMUX = _BV(REFS0) | _BV(MUX3) | _BV(MUX2) | _BV(MUX1);
#endif
  delay(2);                         // Wait for Vref to settle
  ADCSRA |= _BV(ADSC);             // Start conversion
  while (bit_is_set(ADCSRA, ADSC)); // Wait for completion
  result = ADCL;
  result |= ADCH << 8;
  result = 1126400L / result;       // Back-calculate Vcc in mV (1126400 = 1.1 * 1023 * 1000)
  return result;
}

/**
 * Get actual Vcc voltage (with or without correction)
 * @return Vcc in volts
 */
float getVccActual() {
  if (!correct_Vcc) return Vcc_nominal;
  return (readVcc() / 1000.0);      // Convert mV to V
}

/**
 * Read spectrum from AS7341 with automatic gain adjustment
 * Prevents saturation by reducing gain if signal is too strong
 * @param outCh - Array to store 12 channel values (only first 8 are spectral data)
 * @return true if successful, false if read failed
 */
bool safeReadSpectrum(uint16_t outCh[12]) {
  // Initial read with current gain setting
  if (!as.readAllChannels(outCh)) return false;

  // Find maximum value across F1-F8 channels
  uint32_t maxch = 0;
  for (int i = 0; i < 8; ++i) {
    if (outCh[i] > maxch) maxch = outCh[i];
  }

  // If signal is saturating (>60000 counts out of 65535 max), reduce gain
  if (maxch > 60000UL) {
    as.setGain(AS7341_GAIN_1X);     // Reduce to minimum gain
    delay(15);                       // Wait for new gain to take effect
    if (!as.readAllChannels(outCh)) return false;
  }
  
  return true;
}

/**
 * Calculate peak wavelength from spectrum data
 * Uses weighted average method: sum(intensity * wavelength) / sum(intensity)
 * @param ch - Array of 12 channel values
 * @return Peak wavelength in nm, or 0 if signal too weak
 */
float computePeakWavelength(uint16_t ch[12]) {
  double weighted_sum = 0.0;
  double total = 0.0;
  
  // Calculate weighted average across F1-F8 channels
  for (int i = 0; i < 8; ++i) {
    double val = (double) ch[i];
    if (val <= 0.0) continue;        // Skip zero/negative values
    weighted_sum += val * channel_centers_nm[i];
    total += val;
  }
  
  // Return 0 if signal is too weak (< 20 total counts)
  if (total < 20.0) return 0.0;
  
  // Calculate and return weighted average wavelength
  double peak = weighted_sum / total;
  return (float)peak;
}

/**
 * Measure spectrum and send wavelength data to serial
 * Sends two lines:
 * 1. "AS:ch0,ch1,ch2,..." - Raw channel data for debugging
 * 2. "WL,XXX" - Peak wavelength in nm (or 0 if no signal)
 */
void sendSpectrumPeak() {
  uint16_t ch[12];
  
  // Read spectrum with auto-gain adjustment
  if (!safeReadSpectrum(ch)) {
    Serial.println("ERR_SPEC_READ");
    return;
  }

  // Send raw channel data for debugging
  Serial.print("AS:");
  for (int i = 0; i < 12; ++i) {
    Serial.print(ch[i]);
    if (i < 11) Serial.print(',');
  }
  Serial.println();

  // Calculate and send peak wavelength
  float peak_nm = computePeakWavelength(ch);
  if (peak_nm <= 0.0) {
    Serial.println("WL,0");          // No signal detected
  } else {
    Serial.print("WL,");
    Serial.println((int)round(peak_nm)); // Round to nearest nm
  }
}

/**
 * Perform voltage sweep to measure LED I-V curve
 * Sweeps DAC from 0V to maximum, measuring current at each step
 * Outputs CSV format: Vdac(V), I_led(A)
 */
void doSweep() {
  digitalWrite(13, HIGH);            // Turn on onboard LED during sweep
  
  // Send CSV header
  Serial.println("HEADER,Vdac(V),I_led(A)");
  
  // Sweep from 0 to maximum DAC value
  for (int d = 0; d <= DAC_COUNTS; d += sweep_step) {
    // Set DAC output voltage
    dac.setVoltage(d, false);
    delay(settle_ms);                // Wait for voltage to stabilize
    
    // Get current Vcc (supply voltage)
    float vcc_now = getVccActual();
    
    // Calculate DAC output voltage: V_dac = (DAC_value / 4095) * Vcc
    float v_dac = (d / (float)DAC_COUNTS) * vcc_now;
    
    // Measure voltage at LED anode (after resistor)
    float v_led = readAnalogVoltage(measure_led_v, vcc_now);
    
    // Calculate voltage drop across resistor
    float v_res = v_dac - v_led;
    if (v_res < 0.0) v_res = 0.0;    // Clamp to zero if negative (measurement noise)
    
    // Calculate LED current using Ohm's law: I = V / R
    float i_led = v_res / R_series;
    
    // Output data point
    Serial.print(v_dac, 4);          // 4 decimal places
    Serial.print(',');
    Serial.println(i_led, 9);        // 9 decimal places for precision
  }
  
  Serial.println("END");             // Signal end of sweep
  dac.setVoltage(0, false);          // Turn off LED
  digitalWrite(13, LOW);             // Turn off onboard LED
}

/**
 * Turn on LED and measure its spectrum
 * @param dac_value - DAC value to set (0-4095), default is full brightness
 */
void measureLedForSpec(int dac_value = DAC_COUNTS) {
  // Turn on LED at specified brightness
  dac.setVoltage(constrain(dac_value, 0, DAC_COUNTS), false);
  
  // Wait for LED and sensor to stabilize (very important!)
  delay(150);
  
  // Read and send spectrum
  sendSpectrumPeak();
  
  // Turn off LED
  delay(20);
  dac.setVoltage(0, false);
}

// ===== ARDUINO SETUP =====
void setup() {
  Serial.begin(9600);
  Wire.begin();                      // Initialize I2C bus

  // Initialize MCP4725 DAC
  if (!dac.begin(0x60)) {
    // Try alternate address if 0x60 fails
    if (!dac.begin(0x61)) {
      Serial.println("ERR_NO_DAC");
      while (1) { delay(100); }      // Halt if DAC not found
    }
  }

  // Initialize AS7341 spectrometer
  if (!as.begin()) {
    Serial.println("ERR_NO_AS7341");
    while (1) { delay(100); }        // Halt if sensor not found
  }

  // Configure AS7341 sensor
  // Integration time = (ATIME + 1) × (ASTEP + 1) × 2.78µs
  // With ATIME=100, ASTEP=999: ~280ms integration time
  as.setATIME(100);                  // Integration time step count
  as.setASTEP(999);                  // Integration time multiplier
  as.setGain(AS7341_GAIN_8X);        // Gain setting (will auto-adjust if needed)

  // Ensure DAC starts at 0V (LED off)
  dac.setVoltage(0, false);
  
  // Configure onboard LED for status indication
  pinMode(13, OUTPUT);
  
  Serial.println("READY");           // Signal that system is ready
}

// ===== ARDUINO MAIN LOOP =====
void loop() {
  // Check if serial data is available
  if (!Serial.available()) {
    return;
  }
  
  // Read command (terminated by '/')
  serial_command = Serial.readStringUntil('/');
  serial_command.trim();
  if (serial_command.length() == 0) return;

  // ===== COMMAND HANDLERS =====
  
  // Enable Vcc correction
  if (serial_command == "CVcc,t") {
    correct_Vcc = true;
    Serial.println("CVcc,ON");
  } 
  // Disable Vcc correction
  else if (serial_command == "CVcc,f") {
    correct_Vcc = false;
    Serial.println("CVcc,OFF");
  } 
  // Run voltage sweep
  else if (serial_command == "run") {
    doSweep();
  } 
  // Measure LED spectrum and wavelength
  else if (serial_command.startsWith("measure_led")) {
    int v = DAC_COUNTS;              // Default to full brightness
    // Check if custom DAC value was specified
    if (serial_command.indexOf(',') > 0) {
      String s = serial_command.substring(serial_command.indexOf(',') + 1);
      v = s.toInt();
    }
    measureLedForSpec(v);
  } 
  // Unknown command
  else {
    Serial.print("UNK_CMD:");
    Serial.println(serial_command);
  }
}