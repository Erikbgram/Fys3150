import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import imageio as imgio

fig, ax = plt.subplots(figsize=(6,6))

xdata, ydata = [], []
ln, = plt.plot([], [], 'b', label="Earth")
x = []
y = []
z = []

with open("BodyOutput/Earth.txt") as infile:
    infile.readline()
    lines = infile.readlines()
    for line in lines:
        words = line.split(" , ")
        x.append(float(words[0]))
        y.append(float(words[1]))
        z.append(float(words[2]))

x = np.array(x)
y = np.array(y)
z = np.array(z)

def init():
    ax.set_xlim(1.1*min(x), 1.1*max(x))
    ax.set_ylim(1.1*min(y), 1.1*max(y))
    return ln,

def update(frame, frame_x, frame_y, frame_z):
    xdata.append(frame_x[frame])
    ydata.append(frame_y[frame])
    ln.set_data(xdata, ydata)
    #ax.set_xlim(1.1*min(frame_x[:frame+1]), 1.1*max(frame_x[:frame+1]))
    #ax.set_ylim(1.1*min(frame_y[:frame+1]), 1.1*max(frame_y[:frame+1]))
    return ln,

plt.plot(0, 0, "yo", label="Sun")
plt.title("Earth")
plt.xlabel("x")
plt.ylabel("y")
plt.grid()
plt.legend()
plt.tight_layout()

ani = FuncAnimation(fig, update, frames=len(x),
                    fargs=(x,y,z), init_func=init, blit=True, interval=1)
#ani.save("../img/ani.mp4", fps=60)
plt.show()
