import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import imageio as imgio
import sys

dict = {"Sun":"#ffdf22" , "Mercury":"#d3d3d3" , "Venus":"#cde7ff" , "Earth":"#197619" , "Mars":"#cf6e28" , "Jupiter":"#e69500" , "Saturn":"#b9a17e" , "Uranus":"#78d6fb" , "Neptune":"#1414fb" , "Pluto":"#ddb28d"}

filenames = ["forwardEulerbodyOutput_S_E/Sun.txt", "forwardEulerbodyOutput_S_E/Earth.txt"]

"""
x = [[],[]]
y = [[],[]]
z = [[],[]]
total = [[],[]]
data_total = [[],[]]
ln1, = plt.plot([], [], 'b', label="Earth")

for i in range(len(x)):
    with open(filenames[i]) as infile:
        infile.readline()
        lines = infile.readlines()
        for line in lines:
            words = line.split(" , ")
            x[i].append(float(words[0]))
            y[i].append(float(words[1]))
            #z[i].append(float(words[2]))
            #total[i] = [x[i],y[i],z[i]]
            total[i] = [x[i], y[i]]

fig, ax = plt.subplots(figsize=(6,6))

def init():
    ax.set_xlim(1.1*min(x[-1]), 1.1*max(x[-1]))
    ax.set_ylim(1.1*min(y[-1]), 1.1*max(y[-1]))
    return ln1,

def update(frame, frame_total):
    #xdata.append(frame_x[frame])
    #ydata.append(frame_y[frame])
    data_total.append(frame_total)
    ln1.set_data(data_total[2], data_total[1])
    #ax.set_xlim(1.1*min(frame_x[:frame+1]), 1.1*max(frame_x[:frame+1]))
    #ax.set_ylim(1.1*min(frame_y[:frame+1]), 1.1*max(frame_y[:frame+1]))
    return ln1,

plt.plot(0, 0, "yo", label="Sun")
plt.title("Earth")
plt.xlabel("x")
plt.ylabel("y")
plt.grid()
plt.axis("equal")
plt.legend()
plt.tight_layout()

ani = FuncAnimation(fig, update, frames=len(x),
            fargs=(total), init_func=init, blit=True, interval=1)
    #ani.save("../img/ani.gif", writer="pillow", fps=250)
plt.show()
"""


#------------------------------------------------------


import matplotlib.pyplot as plt
from matplotlib import animation
from numpy import random

fig = plt.figure()
#ax1 = plt.axes(xlim=(-108, -104), ylim=(31,34))
ax1 = plt.axes()
line, = ax1.plot([], [], lw=2)
#plt.xlabel('Longitude')
#plt.ylabel('Latitude')

plotlays, plotcols = [2], ["black","red"]
lines = []
for index in range(2):
    lobj = ax1.plot([],[],lw=2,color=plotcols[index])[0]
    lines.append(lobj)

x = [[],[]]
y = [[],[]]
z = [[],[]]
total = [[],[]]
data_total = [[],[]]

def init():

    for i in range(len(x)):
        with open(filenames[i]) as infile:
            infile.readline()
            lines = infile.readlines()
            for line in lines:
                words = line.split(" , ")
                x[i].append(float(words[0]))
                y[i].append(float(words[1]))
                #z[i].append(float(words[2]))
                #total[i] = [x[i],y[i],z[i]]
                total[i] = [x[i], y[i]]

    #print(total[0][0])

    return total

"""
def init():
    for line in lines:
        line.set_data([],[])
    return lines
"""

init()


x1,y1 = [],[]
x2,y2 = [],[]

# fake data
frame_num = 100

def animate(i):

    x = total[0][0, i]
    y = total[1][0, i]
    x1.append(x)
    y1.append(y)

    print(x1)

    x = total[0][1,i]
    y = total[1][1,i]
    x2.append(x)
    y2.append(y)

    xlist = [x1, x2]
    ylist = [y1, y2]

    #for index in range(0,1):
    for lnum,line in enumerate(lines):
        line.set_data(xlist[lnum], ylist[lnum]) # set data for each line separately.

    return lines


# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=frame_num, interval=10, blit=True)


plt.show()



#-----------------------------------------


"""
def eval_planets2(filename1, filename2):

    filenames = [filename1, filename2]

    x = [[],[]]
    y = [[],[]]
    z = [[],[]]
    total = [[],[]]
    data_total = [[],[]]
    ln, = plt.plot([[],[]], [[],[]], 'b', label="Earth")

    for i in range(len(x)):
        with open(filenames[i]) as infile:
            infile.readline()
            lines = infile.readlines()
            for line in lines:
                words = line.split(" , ")
                x[i].append(float(words[0]))
                y[i].append(float(words[1]))
                z[i].append(float(words[2]))
                total[i] = [x[i],y[i],z[i]]

    fig, ax = plt.subplots(figsize=(6,6))

    plt.plot(0, 0, "yo", label="Sun")
    plt.title("Earth")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.grid()
    plt.axis("equal")
    plt.legend()
    plt.tight_layout()

    ani = FuncAnimation(fig, update, frames=len(x),
                        fargs=(total), init_func=init, blit=True, interval=1)
    #ani.save("../img/ani.gif", writer="pillow", fps=250)
    plt.show()
"""


"""
    x = np.array(x)
    y = np.array(y)
    z = np.array(z)

    plt.figure(figsize=(10,6))
    ax = plt.axes()
    # Setting the background color
    ax.set_facecolor("black")
    plt.scatter(x[0], y[0], color = dict["Sun"], label = "Sun")
    plt.scatter(x[1],y[1], color = dict["Earth"], label = "Earth")

    plt.legend()
    plt.title("Two-body solar system")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.grid()
    plt.axis("equal")
    #plt.tight_layout()
    #plt.savefig("../img/earth.png")
    plt.show()
"""

def eval_planets3(filename1, filename2, filename3):

    filenames = [filename1, filename2, filename3]

    x = [[],[],[]]
    y = [[],[],[]]
    z = [[],[],[]]

    for i in range(len(x)):
        with open(filenames[i]) as infile:
            infile.readline()
            lines = infile.readlines()
            for line in lines:
                words = line.split(" , ")
                x[i].append(float(words[0]))
                y[i].append(float(words[1]))
                z[i].append(float(words[2]))

    x = np.array(x)
    y = np.array(y)
    z = np.array(z)

    plt.figure(figsize=(10,6))
    ax = plt.axes()
    # Setting the background color
    ax.set_facecolor("grey")
    plt.scatter(x[0], y[0], color = dict["Sun"], label = "Sun")
    plt.scatter(x[1],y[1], color = dict["Earth"], label = "Earth")
    plt.scatter(x[2],y[2], color = dict["Jupiter"], label = "Jupiter")

    plt.legend()
    plt.title("Three-body solar system")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.grid()
    plt.axis("equal")
    #plt.tight_layout()
    #plt.savefig("../img/earth.png")
    plt.show()

def eval_planets10(filename1, filename2, filename3, filename4, filename5, filename6, filename7, filename8, filename9, filename10):

    filenames = [filename1, filename2, filename3, filename4, filename5, filename6, filename7, filename8, filename9, filename10]

    x = [[],[],[],[],[],[],[],[],[],[]]
    y = [[],[],[],[],[],[],[],[],[],[]]
    z = [[],[],[],[],[],[],[],[],[],[]]
    total = [[],[],[],[],[],[],[],[],[],[]]

    for i in range(len(x)):
        with open(filenames[i]) as infile:
            infile.readline()
            lines = infile.readlines()
            for line in lines:
                words = line.split(" , ")
                x[i].append(float(words[0]))
                y[i].append(float(words[1]))
                z[i].append(float(words[2]))
                total[i] = [[x[i],y[i],z[i]]]

    print(total)

#eval_planets2("forwardEulerbodyOutput_S_E/Sun.txt", "forwardEulerbodyOutput_S_E/Earth.txt")

"""
fig, ax = plt.subplots(figsize=(6,6))

xdata, ydata = [], []
ln, = plt.plot([], [], 'b', label="Earth")
x = []
y = []
z = []
"""
"""
with open("velocityVerletbodyOutput/Earth.txt") as infile:
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
"""


"""
plt.plot(0, 0, "yo", label="Sun")
plt.title("Earth")
plt.xlabel("x")
plt.ylabel("y")
plt.grid()
plt.axis("equal")
plt.legend()
plt.tight_layout()

ani = FuncAnimation(fig, update, frames=len(x),
                    fargs=(x,y,z), init_func=init, blit=True, interval=1)
#ani.save("../img/ani.gif", writer="pillow", fps=250)
plt.show()
"""
