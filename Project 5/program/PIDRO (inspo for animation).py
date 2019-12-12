"""
PID-Regulated Oven
"""
# Imports for general operations
from time import time, localtime, sleep
import logging # Logging module for thermocouple-chip
from os import mkdir

from PID import PID # PID-controller

# Imports for connections
import Adafruit_GPIO # Part of thermocouple configuration
from Adafruit_MAX31856 import MAX31856 as MAX31856 # Library for thermocouple-chip.
from gpiozero import DigitalOutputDevice as DOD # Control for relay

# Imports for plotting
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
from matplotlib.animation import FuncAnimation

# Imports for GUI
import tkinter as tk


init_time = time()
start_time = localtime() # Gives time as [YEAR, MONTH, DAY, HOUR, MINUTE, SECOND, WEEKDAY, YEARDAY, DAYLIGHT SAVINGS]

timestamp = "{:02.0f}{:02.0f}{:.2}-{:02.0f}{:02.0f}".format(start_time[2],
    start_time[1], str(start_time[0])[2:],
    start_time[3], start_time[4]) # Formats time to DDMMYY-hhmm.


try:
    # Create directory for all log-files
    mkdir("logs")
except FileExistsError:
    pass

try:
    # Create directory for MAX31856 log-files
    mkdir("logs/MAX31856")
except FileExistsError:
    pass


with open("logs/general_" + timestamp + ".log", "w") as log:
    pass

logging.basicConfig(filename="logs/MAX31856/" + timestamp + ".log",
    level=logging.DEBUG,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
_logger = logging.getLogger(__name__) # Creates log-file the thermocouple-chip.

# Raspberry Pi hardware SPI configuration.
SPI_PORT   = 0
SPI_DEVICE = 0
sensor = MAX31856(hardware_spi=Adafruit_GPIO.SPI.SpiDev(SPI_PORT, SPI_DEVICE), tc_type=MAX31856.MAX31856_K_TYPE) #remove tc_type for non-K thermocouples

# Relay configuration.
relay = DOD(0) # Creates a general DigitalOutputDevice using GPIO pin 0

# Definitions for plot
fig   = Figure(figsize=(2,2), dpi=100)
ax    = fig.add_subplot(111)
x     = []
y     = []
line, = ax.plot(x,y,"-")

# Animation function
def animate(i):
    temp = sensor.read_temp_c()
    PID_reg(temp)
    elapsed_time = time() - init_time
    x.append(elapsed_time)
    y.append(temp) #float("{:.1f}".format(temp))
    line.set_data(x[:i], y[:i])
    if x[-1]>23:
        del x[0]
        del y[0]
        ax.set_xlim(min(x), max(x))
        ax.set_ylim(min(y)-0.5, max(y)+0.5)
    else:
        ax.set_xlim(min(x), 23)
        ax.set_ylim(min(y)-0.5, max(y)+0.5)

# Setup for PID
P = 1
I = 0
D = 0
pid = PID(P, I, D)
pid.SetPoint = 50
pid.setSampleTime(0.5)

# Choose output using PID
def PID_reg(temp):
    pid.update(temp) # Calculates PID-value
    val = max(min(pid.output, 1),0) # Turns PID-value into something between 0 and 1
    if val<0.5:
        relay.off()
        signal = "Off"
    elif val>=0.5:
        relay.on()
        signal = "On"
    else:
        print("PID-value not recognized. Cutting power.")
        signal = "N/A"
        relay.off()

    with open("logs/general_" + timestamp + ".log", "w") as log:
        log.write("{:.3f}   {:>3}".format(temp, signal))
    print(("Temperature: {:.3f} °C | Signal: {:>3} ").format(temp, signal)) # Prints info

# Placeholder update-function for changing PID
def update(event=None):
    pid.SetPoint = float(sp_ent.get())
    pid.setKp    = float(p_ent.get())
    pid.setKi    = float(i_ent.get())
    pid.setKd    = float(d_ent.get())
    with open("logs/general_" + timestamp + ".log", "w") as log:
        log.write("PID updated: SP: {}, Kp: {}, Ki: {}; Kd: {}".format(pid.SetPoint, pid.setKp, pid.setKi, pid.setKd))
    print("PID updated: Setpoint: {} °C, Kp: {}, Ki: {}; Kd: {}".format(pid.SetPoint, pid.setKp, pid.setKi, pid.setKd))


# Tkinter
window = tk.Tk()
#tk.Grid.rowconfigure(window, 0, weight=1) # Lets us expand Widgets
#tk.Grid.columnconfigure(window, 0, weight=1) # Lets us expand Widgets
window.title("PID-Regulated Oven")
window.configure(background="white")
window.iconbitmap("chip.ico")

# Canvas for matplotlib
canvas = FigureCanvasTkAgg(fig, window)
canvas.get_tk_widget().grid()
canvas.draw()

# Toolbar with grid (not pack)
toolbarFrame = tk.Frame(window)
toolbarFrame.grid(sticky="NSEW")
toolbar = NavigationToolbar2Tk(canvas, toolbarFrame)
toolbar.update()

# Setpoint
tk.Label(window, text="Setpoint", bg="white").grid(row=1, column=2)
sp_ent = tk.Entry(window, width=8, bg="white")
sp_ent.insert(tk.END, "50")
sp_ent.bind("<Return>", update)
sp_ent.grid(row=2, column=2)

# P-value
tk.Label(window, text="P-value", bg="white").grid(row=3, column=1)
p_ent = tk.Entry(window, width=7, bg="white")
p_ent.insert(tk.END, "1")
p_ent.bind("<Return>", update)
p_ent.grid(row=4, column=1)

# I-value
tk.Label(window, text="I-value", bg="white").grid(row=3, column=2)
i_ent = tk.Entry(window, width=7, bg="white")
i_ent.insert(tk.END, "0")
i_ent.bind("<Return>", update)
i_ent.grid(row=4, column=2)

# D-value
tk.Label(window, text="D-value", bg="white").grid(row=3, column=3)
d_ent = tk.Entry(window, textvariable=0, width=7, bg="white")
d_ent.insert(tk.END, "0")
d_ent.bind("<Return>", update)
d_ent.grid(row=4, column=3)



# Starts GUI
ani = FuncAnimation(fig, animate, interval=100)
window.mainloop()
