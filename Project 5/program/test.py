"""

x = [[1,2,3,4,5],[1,2,3,4,5],[1,2,3,4,5]]
y = [[5,4,3,2,1],[5,4,3,2,1],[5,4,3,2,1]]
z = [[1,5,2,4,3], [1,5,2,4,3], [1,5,2,4,3]]

total = [[],[],[]]

for i in range(len(x)):
    total[i] = [x[i],y[i],z[i]]

print(total)
"""

#----------------------- Animasjon 1 ------------------------

"""
import matplotlib.pyplot as plt
from matplotlib import animation
from numpy import random

fig = plt.figure()
ax1 = plt.axes(xlim=(-108, -104), ylim=(31,34))
line, = ax1.plot([], [], lw=2)
plt.xlabel('Longitude')
plt.ylabel('Latitude')

plotlays, plotcols = [2], ["black","red"]
lines = []
for index in range(2):
    lobj = ax1.plot([],[],lw=2,color=plotcols[index])[0]
    lines.append(lobj)


def init():
    for line in lines:
        line.set_data([],[])
    return lines

x1,y1 = [],[]
x2,y2 = [],[]

# fake data
frame_num = 100
gps_data = [-104 - (4 * random.rand(2, frame_num)), 31 + (3 * random.rand(2, frame_num))]


def animate(i):

    x = gps_data[0][0, i]
    y = gps_data[1][0, i]
    x1.append(x)
    y1.append(y)

    x = gps_data[0][1,i]
    y = gps_data[1][1,i]
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
"""

#----------------------- Animasjon 2 -----------------------

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
plt.style.use('dark_background')

fig = plt.figure()
ax = plt.axes(xlim=(-50, 50), ylim=(-50, 50))
line, = ax.plot([], [], lw=2)

# initialization function
def init():
	# creating an empty plot/frame
	line.set_data([], [])
	return line,

# lists to store x and y axis points
xdata, ydata = [], []

# animation function
def animate(i):
	# t is a parameter
	t = 0.1*i

	# x, y values to be plotted
	x = t*np.sin(t)
	y = t*np.cos(t)

	# appending new points to x, y axes points list
	xdata.append(x)
	ydata.append(y)
	line.set_data(xdata, ydata)
	return line,

# setting a title for the plot
plt.title('Creating a growing coil with matplotlib!')
# hiding the axis details
plt.axis('off')

# call the animator
anim = animation.FuncAnimation(fig, animate, init_func=init,
							frames=500, interval=20, blit=True)

plt.show()

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

#filenames = ["forwardEulerbodyOutput_S_E/Sun.txt", "forwardEulerbodyOutput_S_E/Earth.txt"]

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


#-------------------------------- Animasjon 1 ----------------------------

"""
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
lines, = []
for index in range(2):
    lobj = ax1.plot([],[],lw=2,color=plotcols[index])[0]
    lines.append(lobj)


def init():
    for line in lines:
        line.set_data([],[])
    return lines,

x = [[],[]]
y = [[],[]]
z = [[],[]]
total = [[],[]]
data_total = [[],[]]

x1,y1 = [],[]
x2,y2 = [],[]

# fake data
frame_num = 100

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


def animate(i):

    X1 = total[0][0,i]
    Y1 = total[1][0,i]
    x1.append(X1)
    y1.append(Y1)

    X2 = total[0][1,i]
    Y2 = total[1][1,i]
    x2.append(X2)
    y2.append(Y2)


    xlist = [x1, x2]
    ylist = [y1, y2]

    #for index in range(0,1):
    for lnum,line in enumerate(lines):
        line.set_data(xlist[lnum], ylist[lnum]) # set data for each line separately.

    return lines,


# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=frame_num, interval=10, blit=True)


plt.show()
"""


#------------------------------------ Animasjon 2 -----------------------------------

"""
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
#plt.style.use('dark_background')

fig = plt.figure()
ax = plt.axes(xlim=(-50, 50), ylim=(-50, 50))
line, = ax.plot([], [], lw=2)

# initialization function
def init():
	# creating an empty plot/frame
	line.set_data([], [])
	return line,

# lists to store x and y axis points
xdata, ydata = [], []

# animation function
def animate(i):

    x = [[],[]]
    y = [[],[]]
    z = [[],[]]
    total = [[],[]]
    data_total = [[],[]]

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
                #total[i] = [x[i], y[i]]


    #x = total[1][1][0]
    #y = total[1][1][1]

    # appending new points to x, y axes points list
    #xdata.append(x)
    #ydata.append(y)
    #line.set_data(xdata, ydata)
    line.set_data(x[0], y[0])
    return line,

# setting a title for the plot
plt.title('Creating a growing coil with matplotlib!')
# hiding the axis details
plt.axis('off')

# call the animator
anim = animation.FuncAnimation(fig, animate, init_func=init,
							frames=500, interval=20, blit=True)

plt.show()

"""

#--------------------------- Animasjon 3 ----------------------------------
# DENNE FUNGERER!!!!! - den som blir brukt i animasjonsfilene
