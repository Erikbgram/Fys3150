import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import mpl_toolkits.mplot3d.axes3d as p3

fig = plt.figure()
ax = p3.Axes3D(fig)
xdata, ydata, zdata = [], [], []
ln, = plt.plot([], [], 'b', label="Earth")
t = []
x = []
y = []
z = []

with open("output.txt") as infile:
    infile.readline()
    lines = infile.readlines()
    for line in lines:
        words = line.split(", ")
        t.append(float(words[0]))
        x.append(float(words[1]))
        y.append(float(words[2]))
        z.append(float(words[3]))

t = np.array(t)
x = np.array(x)
y = np.array(y)
z = np.array(z)

data = 

def init():
    ax.set_xlim(1.1*min(x), 1.1*max(x))
    ax.set_ylim(1.1*min(y), 1.1*max(y))
    ax.set_zlim(1.1*min(z), 1.1*max(z))
    return ln,

def update(frame, frame_t, frame_x, frame_y, frame_z):
    xdata.append(frame_x[frame])
    ydata.append(frame_y[frame])
    zdata.append(frame_z[frame])
    ln.set_data(xdata, ydata, zdata)
    ln.set_3d_properties(data[2, :frame])
    #ax.set_xlim(1.1*min(frame_x[:frame+1]), 1.1*max(frame_x[:frame+1]))
    #ax.set_ylim(1.1*min(frame_y[:frame+1]), 1.1*max(frame_y[:frame+1]))
    return ln,

#plt.plot(0, 0, 0, "yo", label="Sun")
plt.title("Earth")
plt.xlabel("x")
plt.ylabel("y")
plt.ylabel("z")
plt.grid()
plt.legend()

ani = FuncAnimation(fig, update, frames=len(t),
                    fargs=(t,x,y,z), init_func=init, blit=True, interval=1)
#ani.save("../img/ani.mp4", fps=60)
plt.show()
