"""
LED Planck Constant Measurement System - Python GUI

This program:
1. Controls Arduino to sweep LED voltage and measure I-V curve
2. Estimates LED threshold voltage (bandgap voltage) from I-V data
3. Measures LED peak wavelength using spectrometer
4. Saves wavelength and threshold voltage pairs for multiple LEDs
5. Calculates Planck's constant from the relationship: eV_th = hc/λ

Theory:
- When LED turns on at threshold voltage V_th, the photon energy equals bandgap: eV_th = E_g
- Photon energy also equals: E = hc/λ
- Therefore: eV_th = hc/λ
- Rearranging: V_th = (hc/e) * (1/λ)
- This is a linear relationship where slope = hc/e
- From slope, we can calculate: h = slope * e / c

Where:
- h = Planck's constant (6.626 × 10^-34 J·s)
- c = speed of light (3 × 10^8 m/s)
- e = elementary charge (1.602 × 10^-19 C)
- λ = wavelength (m)
- V_th = threshold voltage (V)
"""

import serial
import threading
import time
import numpy as np
import tkinter as tk
from tkinter import ttk, messagebox
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from scipy.stats import linregress
import csv
import os

# ===== CONFIGURATION =====
BAUD = 9600                          # Serial communication baud rate
CSV_FILE = "led_planck_data.csv"     # File to store LED measurements

class LEDGUI:
    """Main GUI application for LED Planck constant measurement"""
    
    def __init__(self, root):
        """
        Initialize the GUI application
        
        Args:
            root: Tkinter root window
        """
        self.root = root
        root.title("LED Planck Constant Lab")
        
        # Serial connection
        self.ser = None
        self.port_var = tk.StringVar(value="COM3")  # Default serial port
        self.is_sweeping = False                     # Flag for sweep in progress
        
        # Data storage
        self.vdac = []              # DAC voltages from sweep
        self.i_led = []             # LED currents from sweep
        self.wavelengths = []       # Measured wavelengths (nm)
        self.thresholds = []        # Calculated threshold voltages (V)
        
        # ===== BUILD USER INTERFACE =====
        
        # Control panel frame
        frame = ttk.Frame(root, padding=8)
        frame.pack(fill="x")
        
        # Serial port connection
        ttk.Label(frame, text="Serial Port:").grid(row=0, column=0, sticky="w")
        ttk.Entry(frame, textvariable=self.port_var, width=12).grid(row=0, column=1)
        ttk.Button(frame, text="Connect", command=self.connect).grid(row=0, column=2, padx=6)
        
        # Action buttons
        ttk.Button(frame, text="1. Run Sweep", command=self.start_sweep).grid(row=1, column=0, pady=6, sticky="ew")
        ttk.Button(frame, text="2. Measure Wavelength", command=self.measure_spec).grid(row=1, column=1, sticky="ew")
        ttk.Button(frame, text="3. Compute Planck's h", command=self.compute_planck).grid(row=1, column=2, sticky="ew")
        
        # Clear data button
        ttk.Button(frame, text="Clear All Data", command=self.clear_data).grid(row=2, column=0, columnspan=3, pady=6)
        
        # Status label
        self.status_label = ttk.Label(frame, text="Status: Not connected", foreground="red")
        self.status_label.grid(row=3, column=0, columnspan=3, pady=6)
        
        # ===== MATPLOTLIB PLOT =====
        
        # Create figure and axis for I-V curve
        self.fig, self.ax = plt.subplots(figsize=(8, 5))
        self.ax.set_xlabel("V_dac (V)", fontsize=12)
        self.ax.set_ylabel("I_led (A)", fontsize=12)
        self.ax.set_title("LED I-V Characteristic Curve", fontsize=14)
        self.ax.grid(True, alpha=0.3)
        
        # Create empty plot line
        self.line, = self.ax.plot([], [], '-o', markersize=3, linewidth=1.5)
        
        # Embed matplotlib figure in tkinter window
        self.canvas = FigureCanvasTkAgg(self.fig, master=root)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
    def connect(self):
        """
        Connect to Arduino via serial port
        Opens serial connection and waits for Arduino to initialize
        """
        port = self.port_var.get()
        try:
            # Close existing connection if open
            if self.ser and self.ser.is_open:
                self.ser.close()
            
            # Open new serial connection
            self.ser = serial.Serial(port, BAUD, timeout=1)
            time.sleep(2)  # Wait for Arduino to reset and initialize
            
            # Update status
            self.status_label.config(text=f"Status: Connected to {port}", foreground="green")
            messagebox.showinfo("Serial", f"Connected to {port}")
            
        except Exception as e:
            self.status_label.config(text="Status: Connection failed", foreground="red")
            messagebox.showerror("Serial", f"Failed to open {port}: {e}")
    
    def start_sweep(self):
        """
        Start voltage sweep measurement
        Runs sweep in separate thread to avoid freezing GUI
        """
        # Check if serial is connected
        if not self.ser or not self.ser.is_open:
            messagebox.showerror("Error", "Serial not connected. Click 'Connect' first.")
            return
        
        # Check if sweep already running
        if self.is_sweeping:
            messagebox.showwarning("Warning", "Sweep already in progress")
            return
        
        # Clear previous data
        self.vdac.clear()
        self.i_led.clear()
        
        # Start sweep in background thread
        t = threading.Thread(target=self._sweep_thread, daemon=True)
        t.start()
        
        self.status_label.config(text="Status: Running voltage sweep...", foreground="blue")
    
    def _sweep_thread(self):
        """
        Background thread that performs voltage sweep
        Reads data from Arduino and updates plot in real-time
        """
        self.is_sweeping = True
        
        try:
            # Clear serial buffer
            self.ser.reset_input_buffer()
            
            # Send sweep command to Arduino
            self.ser.write(b"run/")
            
            # Read data until "END" received
            while True:
                line = self.ser.readline().decode(errors='ignore').strip()
                
                if not line:
                    continue  # Skip empty lines
                
                if line.startswith("HEADER"):
                    continue  # Skip header line
                
                if line == "END":
                    break  # Sweep complete
                
                # Parse voltage and current data
                try:
                    parts = line.split(',')
                    if len(parts) >= 2:
                        v = float(parts[0])  # Voltage
                        i = float(parts[1])  # Current
                        
                        # Store data
                        self.vdac.append(v)
                        self.i_led.append(i)
                        
                        # Update plot (every 10 points for performance)
                        if len(self.vdac) % 10 == 0:
                            self.line.set_data(self.vdac, self.i_led)
                            self.ax.relim()
                            self.ax.autoscale_view()
                            self.canvas.draw_idle()
                
                except (ValueError, IndexError):
                    # Ignore malformed lines
                    pass
            
            # Final plot update
            self.line.set_data(self.vdac, self.i_led)
            self.ax.relim()
            self.ax.autoscale_view()
            self.canvas.draw()
            
            # Estimate threshold voltage from I-V data
            thresh = self.estimate_threshold()
            
            if thresh is not None:
                self.thresholds.append(thresh)
                self.status_label.config(
                    text=f"Status: Sweep complete. V_threshold = {thresh:.4f} V", 
                    foreground="green"
                )
                messagebox.showinfo("Threshold Voltage", 
                    f"Estimated threshold voltage:\n\nV_th = {thresh:.4f} V\n\n"
                    f"This is the LED's bandgap voltage.\n"
                    f"Now click 'Measure Wavelength' to measure the LED color."
                )
            else:
                self.status_label.config(text="Status: Could not estimate threshold", foreground="orange")
                messagebox.showwarning("Threshold Voltage", 
                    "Could not reliably estimate threshold voltage.\n"
                    "The I-V curve may be too noisy or the LED is not conducting.\n"
                    "Check your wiring and try again."
                )
        
        except Exception as e:
            messagebox.showerror("Error", f"Sweep failed: {e}")
            self.status_label.config(text="Status: Sweep failed", foreground="red")
        
        finally:
            self.is_sweeping = False
    
    def estimate_threshold(self):
        """
        Estimate LED threshold voltage from I-V curve
        
        Method: Linear extrapolation
        1. Find the region where current is rising (above noise floor)
        2. Fit a line to this linear region
        3. Find x-intercept (where I = 0) - this is V_threshold
        
        Returns:
            float: Threshold voltage in volts, or None if estimation fails
        """
        v = np.array(self.vdac)
        i = np.array(self.i_led)
        
        # Need minimum data points
        if len(v) < 10:
            return None
        
        # Estimate noise floor from first 10% of data
        noise = np.median(i[:max(3, len(i)//10)])
        
        # Find points significantly above noise (linear region)
        mask = i > max(noise * 4, 1e-6)  # 4x noise or 1µA minimum
        
        # If not enough points, try looser threshold
        if mask.sum() < 6:
            mask = i > max(noise * 2, 1e-7)
        
        # If still not enough, use second half of data
        if mask.sum() < 6:
            mask = np.arange(len(i)) > (len(i) // 2)
        
        # Get masked data points
        xv = v[mask]
        yv = i[mask]
        
        if len(xv) < 6:
            return None
        
        # Fit line to linear region: I = slope * V + intercept
        slope, intercept, r, p, se = linregress(xv, yv)
        
        # Check if slope is reasonable
        if abs(slope) < 1e-12:
            return None
        
        # Calculate x-intercept where I = 0
        # 0 = slope * V_th + intercept
        # V_th = -intercept / slope
        V_th = -intercept / slope
        
        # Sanity check: threshold should be positive and reasonable
        if V_th < 0 or V_th > max(v) * 1.2:
            return None
        
        return float(V_th)
    
    def measure_spec(self):
        """
        Measure LED wavelength using spectrometer
        Turns on LED and reads spectrum from AS7341 sensor
        """
        # Check if serial is connected
        if not self.ser or not self.ser.is_open:
            messagebox.showerror("Error", "Serial not connected")
            return
        
        # Check if we have a threshold measurement
        if len(self.thresholds) == 0:
            response = messagebox.askywarning("Warning", 
                "No threshold voltage measured yet.\n"
                "You should run a voltage sweep first.\n\n"
                "Continue anyway?"
            )
            if not response:
                return
        
        self.status_label.config(text="Status: Measuring wavelength...", foreground="blue")
        
        try:
            # Clear serial buffer
            self.ser.reset_input_buffer()
            
            # Send measure command to Arduino
            self.ser.write(b"measure_led/")
            
            peak_nm = None
            
            # Read until we get wavelength or timeout
            t0 = time.time()
            while time.time() - t0 < 3.0:  # 3 second timeout
                line = self.ser.readline().decode(errors='ignore').strip()
                
                if not line:
                    continue
                
                # Skip debug output
                if line.startswith("AS:"):
                    continue
                
                # Look for wavelength line
                if line.startswith("WL,"):
                    try:
                        peak_nm = int(line.split(",")[1])
                        break
                    except (ValueError, IndexError):
                        pass
            
            # Check if we got a wavelength
            if peak_nm is None or peak_nm == 0:
                self.status_label.config(text="Status: No wavelength detected", foreground="orange")
                messagebox.showwarning("Spectrometer", 
                    "No wavelength detected.\n\n"
                    "Possible issues:\n"
                    "- LED is too dim or not on\n"
                    "- LED is not pointing at spectrometer\n"
                    "- Spectrometer is not working\n\n"
                    "Try positioning LED closer to spectrometer."
                )
                return
            
            # Store wavelength
            self.wavelengths.append(peak_nm)
            
            # If we have a corresponding threshold, save to CSV
            if len(self.thresholds) >= 1:
                V_th = self.thresholds[-1]
                self.append_csv(peak_nm, V_th)
                
                self.status_label.config(
                    text=f"Status: Saved λ={peak_nm}nm, V_th={V_th:.4f}V", 
                    foreground="green"
                )
                
                messagebox.showinfo("Wavelength Measured", 
                    f"Peak wavelength: {peak_nm} nm\n"
                    f"Threshold voltage: {V_th:.4f} V\n\n"
                    f"Data saved to {CSV_FILE}\n\n"
                    f"You can now:\n"
                    f"1. Measure another LED (repeat steps 1-2)\n"
                    f"2. Calculate Planck's constant (click button 3)"
                )
            else:
                self.status_label.config(text=f"Status: Wavelength={peak_nm}nm (no V_th)", foreground="orange")
                messagebox.showinfo("Wavelength Measured", 
                    f"Peak wavelength: {peak_nm} nm\n\n"
                    f"No threshold voltage available.\n"
                    f"Run a voltage sweep first to get complete data."
                )
        
        except Exception as e:
            messagebox.showerror("Error", f"Wavelength measurement failed: {e}")
            self.status_label.config(text="Status: Measurement failed", foreground="red")
    
    def append_csv(self, wavelength_nm, v_threshold):
        """
        Append LED measurement to CSV file
        
        Args:
            wavelength_nm: Peak wavelength in nanometers
            v_threshold: Threshold voltage in volts
        """
        # Check if file exists (to determine if we need header)
        header_needed = not os.path.exists(CSV_FILE)
        
        # Append data to CSV
        with open(CSV_FILE, 'a', newline='') as f:
            writer = csv.writer(f)
            
            # Write header if new file
            if header_needed:
                writer.writerow(["wavelength_nm", "V_th_V"])
            
            # Write data
            writer.writerow([wavelength_nm, v_threshold])
    
    def compute_planck(self):
        """
        Calculate Planck's constant from collected LED data
        
        Theory:
        - eV_th = hc/λ  (photon energy at threshold)
        - Rearranging: V_th = (hc/e) × (1/λ)
        - This is linear: V_th = slope × (1/λ) + intercept
        - slope = hc/e
        - Therefore: h = slope × e / c
        """
        # Check if data file exists
        if not os.path.exists(CSV_FILE):
            messagebox.showerror("Error", 
                f"No saved LED data found.\n\n"
                f"File '{CSV_FILE}' does not exist.\n"
                f"Please measure at least 3 different LEDs first."
            )
            return
        
        try:
            # Load data from CSV
            data = np.loadtxt(CSV_FILE, delimiter=',', skiprows=1)
            
            # Handle single data point case
            if data.ndim == 1 and len(data) == 2:
                data = data.reshape((1, 2))
            
            # Check if we have enough data points
            if len(data) < 3:
                messagebox.showwarning("Insufficient Data", 
                    f"Only {len(data)} LED(s) measured.\n\n"
                    f"You need at least 3 different colored LEDs\n"
                    f"for a reliable Planck's constant calculation.\n\n"
                    f"Measure more LEDs and try again."
                )
                return
            
            # Extract wavelength and threshold voltage
            lam_nm = data[:, 0]     # Wavelength in nm
            Vth = data[:, 1]        # Threshold voltage in V
            
            # Convert wavelength to meters
            lam_m = lam_nm * 1e-9
            
            # Calculate 1/λ for linear fit
            inv_lam = 1.0 / lam_m
            
            # Linear regression: V_th = slope × (1/λ) + intercept
            slope, intercept, r, p, se = linregress(inv_lam, Vth)
            
            # Calculate Planck's constant
            # slope = hc/e
            # h = slope × e / c
            e = 1.602176634e-19     # Elementary charge (C)
            c = 2.99792458e8        # Speed of light (m/s)
            h_est = slope * e / c   # Estimated Planck's constant (J·s)
            
            # Known value for comparison
            h_actual = 6.62607015e-34  # J·s
            error_percent = abs((h_est - h_actual) / h_actual) * 100
            
            # Display results
            msg = (
                f"Planck's Constant Calculation Results\n"
                f"{'='*45}\n\n"
                f"Number of LEDs measured: {len(data)}\n\n"
                f"Linear fit: V_th = m × (1/λ) + b\n"
                f"  Slope (m) = {slope:.4e} V·m\n"
                f"  Intercept (b) = {intercept:.4f} V\n"
                f"  R² = {r**2:.4f}\n\n"
                f"Estimated Planck's constant:\n"
                f"  h = {h_est:.4e} J·s\n\n"
                f"Accepted value:\n"
                f"  h = {h_actual:.4e} J·s\n\n"
                f"Error: {error_percent:.1f}%\n\n"
                f"{'='*45}\n"
                f"Interpretation:\n"
                f"  R² > 0.95: Excellent fit\n"
                f"  Error < 10%: Good result\n"
                f"  Error < 20%: Acceptable for lab"
            )
            
            self.status_label.config(
                text=f"Status: h = {h_est:.3e} J·s (error: {error_percent:.1f}%)", 
                foreground="green"
            )
            
            messagebox.showinfo("Planck's Constant Result", msg)
            
            # Create plot of fit
            self._plot_planck_fit(inv_lam, Vth, slope, intercept, lam_nm, h_est, error_percent)
        
        except Exception as e:
            messagebox.showerror("Error", f"Calculation failed: {e}")
            self.status_label.config(text="Status: Calculation failed", foreground="red")
    
    def _plot_planck_fit(self, inv_lam, Vth, slope, intercept, lam_nm, h_est, error_percent):
        """
        Create plot showing linear fit for Planck's constant calculation
        
        Args:
            inv_lam: 1/wavelength values (1/m)
            Vth: Threshold voltages (V)
            slope: Fit slope (V·m)
            intercept: Fit intercept (V)
            lam_nm: Wavelengths in nm (for color coding)
            h_est: Estimated Planck's constant
            error_percent: Error percentage
        """
        # Create new figure
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Color map for LED wavelengths
        # Red LEDs = red, yellow = yellow, green = green, blue = blue
        colors = []
        for wl in lam_nm:
            if wl < 490:
                colors.append('blue')
            elif wl < 530:
                colors.append('cyan')
            elif wl < 570:
                colors.append('green')
            elif wl < 600:
                colors.append('yellow')
            else:
                colors.append('red')
        
        # Scatter plot of data points
        for i in range(len(inv_lam)):
            ax.scatter(inv_lam[i], Vth[i], s=100, c=colors[i], 
                      edgecolor='black', linewidth=2, zorder=3,
                      label=f'{int(lam_nm[i])} nm')
        
        # Plot fit line
        x_fit = np.linspace(min(inv_lam) * 0.9, max(inv_lam) * 1.1, 100)
        y_fit = slope * x_fit + intercept
        ax.plot(x_fit, y_fit, 'k--', linewidth=2, label='Linear fit', zorder=2)
        
        # Labels and title
        ax.set_xlabel('1 / λ  (m⁻¹)', fontsize=14)
        ax.set_ylabel('Threshold Voltage V_th (V)', fontsize=14)
        ax.set_title(
            f"Planck's Constant from LED Data\n"
            f"h = {h_est:.3e} J·s  (error: {error_percent:.1f}%)",
            fontsize=16
        )
        ax.grid(True, alpha=0.3)
        ax.legend(loc='best', fontsize=10)
        
        # Add equation to plot
        textstr = f'V_th = ({slope:.2e})×(1/λ) + {intercept:.3f}\nR² = {np.corrcoef(inv_lam, Vth)[0,1]**2:.4f}'
        ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=12,
                verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        plt.tight_layout()
        plt.show()
    
    def clear_data(self):
        """Clear all collected data and CSV file"""
        response = messagebox.askyesno("Clear Data", 
            f"This will delete:\n"
            f"- Current session data\n"
            f"- {CSV_FILE}\n\n"
            f"Are you sure?"
        )
        
        if response:
            # Clear in-memory data
            self.vdac.clear()
            self.i_led.clear()
            self.wavelengths.clear()
            self.thresholds.clear()
            
            # Clear plot
            self.line.set_data([], [])
            self.ax.relim()
            self.ax.autoscale_view()
            self.canvas.draw()
            
            # Delete CSV file
            if os.path.exists(CSV_FILE):
                os.remove(CSV_FILE)
            
            self.status_label.config(text="Status: All data cleared", foreground="blue")
            messagebox.showinfo("Data Cleared", "All data has been cleared.")

def main():
    """Main entry point for the application"""
    root = tk.Tk()
    app = LEDGUI(root)
    root.mainloop()

if __name__ == "__main__":
    main()