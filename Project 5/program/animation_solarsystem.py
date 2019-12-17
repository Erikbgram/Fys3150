#Sist endret: 16.12.2019 21:00 by Alexandra Jahr Kolstad

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
#import imageio as imgio
#import sys

dict = {"Sun":"#ffdf22" , "Mercury":"#d3d3d3" , "Venus":"#c04e01" , "Earth":"#197619" , "Mars":"#cf6e28" , "Jupiter":"#e69500" , "Saturn":"#b9a17e" , "Uranus":"#78d6fb" , "Neptune":"#1414fb" , "Pluto":"#ddb28d"}

def eval_planets2(filename1, filename2):

    filenames = [filename1, filename2]

    x = [[],[]]
    y = [[],[]]
    z = [[],[]]

    for i in range(len(x)):
        with open(filenames[i]) as infile:
            infile.readline()
            lines = infile.readlines()
            for line in lines:
                words = line.split(" , ")
                x[i].append(float(words[0]))
                y[i].append(float(words[1]))

    fig, ax = plt.subplots()
    line1, = plt.plot(x[0], y[0], color = dict["Sun"], label = "Sun")
    line2, = plt.plot(x[1], y[1], color = dict["Earth"], label = "Earth")

    ani = animation.FuncAnimation(fig, update_planets2, frames = 500, fargs=[x, y, line1, line2],
                      interval=10, blit=True)

    ax.set_facecolor("black")           # Setting the background color
    plt.legend()
    plt.title("Two-body solar system with " + filename1[:-22])
    plt.xlabel("x")
    plt.ylabel("y")
    plt.grid()
    plt.axis("equal")
    plt.tight_layout()
    #plt.savefig("../img/earth.png")
    plt.show()

def update_planets2(num, x, y, line1, line2):
    line1.set_data(x[0][:num], y[0][:num])
    line2.set_data(x[1][:num], y[1][:num])
    return [line1,line2]

def eval_planets3(filename1, filename2, filename3):

    filenames = [filename1, filename2, filename3]

    x = [[],[],[]]
    y = [[],[],[]]

    for i in range(len(x)):
        with open(filenames[i]) as infile:
            infile.readline()
            lines = infile.readlines()
            for line in lines:
                words = line.split(" , ")
                x[i].append(float(words[0]))
                y[i].append(float(words[1]))


    fig, ax = plt.subplots()
    line1, = plt.plot(x[0], y[0], color = dict["Sun"], label = "Sun")
    line2, = plt.plot(x[1], y[1], color = dict["Earth"], label = "Earth")
    line3, = plt.plot(x[2], y[2], color = dict["Jupiter"], label = "Jupiter")

    ani = animation.FuncAnimation(fig, update_planets3, frames = 500, fargs=[x, y, line1, line2, line3],
                      interval=10, blit=True)


    ax.set_facecolor("grey")            # Setting the background color
    plt.legend()
    plt.title("Three-body solar system with " + filename1[:-24])
    plt.xlabel("x")
    plt.ylabel("y")
    plt.grid()
    plt.axis("equal")
    plt.tight_layout()
    #plt.savefig("../img/earth.png")
    plt.show()

def update_planets3(num, x, y, line1, line2, line3):
    line1.set_data(x[0][:num], y[0][:num])
    line2.set_data(x[1][:num], y[1][:num])
    line3.set_data(x[2][:num], y[2][:num])
    return [line1,line2,line3]

def eval_planets10(filename1, filename2, filename3, filename4, filename5, filename6, filename7, filename8, filename9, filename10):

    filenames = [filename1, filename2, filename3, filename4, filename5, filename6, filename7, filename8, filename9, filename10]

    x = [[],[],[],[],[],[],[],[],[],[]]
    y = [[],[],[],[],[],[],[],[],[],[]]

    for i in range(len(x)):
        with open(filenames[i]) as infile:
            infile.readline()
            lines = infile.readlines()
            for line in lines:
                words = line.split(" , ")
                x[i].append(float(words[0]))
                y[i].append(float(words[1]))

    fig, ax = plt.subplots()
    line1, = plt.plot(x[0], y[0], color = dict["Sun"], label = "Sun")
    line2, = plt.plot(x[1], y[1], color = dict["Mercury"], label = "Mercury")
    line3, = plt.plot(x[2], y[2], color = dict["Venus"], label = "Venus")
    line4, = plt.plot(x[3], y[3], color = dict["Earth"], label = "Earth")
    line5, = plt.plot(x[4], y[4], color = dict["Mars"], label = "Mars")
    line6, = plt.plot(x[5], y[5], color = dict["Jupiter"], label = "Jupiter")
    line7, = plt.plot(x[6], y[6], color = dict["Saturn"], label = "Saturn")
    line8, = plt.plot(x[7], y[7], color = dict["Uranus"], label = "Uranus")
    line9, = plt.plot(x[8], y[8], color = dict["Neptune"], label = "Neptune")
    line10, = plt.plot(x[9], y[9], color = dict["Pluto"], label = "Pluto")

    ani = animation.FuncAnimation(fig, update_planets10, frames = 500, fargs=[x, y, line1, line2, line3, line4, line5, line6, line7, line8, line9, line10],
                      interval=10, blit=True)

    ax.set_facecolor("black")            # Setting the background color
    plt.legend(loc = "center left")
    #ax = plt.axes(xlim=(-150, 150))
    #plt.axis([-150,150,-150,150])
    #plt.legend()
    plt.title("Ten-body solar system with " + filename1[:-18])
    plt.xlabel("x")
    plt.ylabel("y")

    plt.grid()
    plt.axis("equal")
    plt.tight_layout()
    #plt.savefig("../img/earth.png")
    plt.show()

def update_planets10(num, x, y, line1, line2, line3, line4, line5, line6, line7, line8, line9, line10):
    line1.set_data(x[0][:num], y[0][:num])
    line2.set_data(x[1][:num], y[1][:num])
    line3.set_data(x[2][:num], y[2][:num])
    line4.set_data(x[3][:num], y[3][:num])
    line5.set_data(x[4][:num], y[4][:num])
    line6.set_data(x[5][:num], y[5][:num])
    line7.set_data(x[6][:num], y[6][:num])
    line8.set_data(x[7][:num], y[7][:num])
    line9.set_data(x[8][:num], y[8][:num])
    line10.set_data(x[9][:num], y[9][:num])
    return [line1,line2,line3,line4,line5,line6,line7,line8,line9,line10]


#-------------------------------------------

"""
#ny

#Forward Euler output

eval_planets2("ForwardEulerbodyOutput_S_E/Sun.txt", "ForwardEulerbodyOutput_S_E/Earth.txt")

#eval_planets3("ForwardEulerbodyOutput_S_E_J/Sun.txt", "ForwardEulerbodyOutput_S_E_J/Earth.txt", "ForwardEulerbodyOutput_S_E_J/Jupiter.txt")

eval_planets10("ForwardEulerbodyOutput_10body/Sun.txt", "ForwardEulerbodyOutput_10body/Mercury.txt", "ForwardEulerbodyOutput_10body/Venus.txt", "ForwardEulerbodyOutput_10body/Earth.txt", "ForwardEulerbodyOutput_10body/Mars.txt", "ForwardEulerbodyOutput_10body/Jupiter.txt", "ForwardEulerbodyOutput_10body/Saturn.txt", "ForwardEulerbodyOutput_10body/Uranus.txt", "ForwardEulerbodyOutput_10body/Neptune.txt", "ForwardEulerbodyOutput_10body/Pluto.txt")

#Velocity verlet output

#eval_planets2("VelocityVerletbodyOutput_S_E/Sun.txt", "VelocityVerletbodyOutput_S_E/Earth.txt")

#eval_planets3(""VelocityVerletbodyOutput_S_E_J/Sun.txt", "VelocityVerletbodyOutput_S_E_J/Earth.txt", "VelocityVerletbodyOutput_S_E_J/Jupiter.txt")

eval_planets10("VelocityVerletbodyOutput_10body/Sun.txt", "VelocityVerletbodyOutput_10body/Mercury.txt", "VelocityVerletbodyOutput_10body/Venus.txt", "VelocityVerletbodyOutput_10body/Earth.txt", "VelocityVerletbodyOutput_10body/Mars.txt", "VelocityVerletbodyOutput_10body/Jupiter.txt", "VelocityVerletbodyOutput_10body/Saturn.txt", "VelocityVerletbodyOutput_10body/Uranus.txt", "VelocityVerletbodyOutput_10body/Neptune.txt", "VelocityVerletbodyOutput_10body/Pluto.txt")
"""


#gammel

#Forward Euler output

eval_planets2("ForwardEulerbodyOutput_S_E/Sun.txt", "ForwardEulerbodyOutput_S_E/Earth.txt")

#eval_planets3("ForwardEulerbodyOutput_S_E_J/Sun.txt", "ForwardEulerbodyOutput_S_E_J/Earth.txt", "ForwardEulerbodyOutput_S_E_J/Jupiter.txt")

eval_planets10("ForwardEulerbodyOutput_10body/Sun.txt", "ForwardEulerbodyOutput_10body/Mercury.txt", "ForwardEulerbodyOutput_10body/Venus.txt", "ForwardEulerbodyOutput_10body/Earth.txt", "ForwardEulerbodyOutput_10body/Mars.txt", "ForwardEulerbodyOutput_10body/Jupiter.txt", "ForwardEulerbodyOutput_10body/Saturn.txt", "ForwardEulerbodyOutput_10body/Uranus.txt", "ForwardEulerbodyOutput_10body/Neptune.txt", "ForwardEulerbodyOutput_10body/Pluto.txt")

#Velocity verlet output

#eval_planets2("VelocityVerletbodyOutput_S_E/Sun.txt", "VelocityVerletbodyOutput_S_E/Earth.txt")

#eval_planets3(""VelocityVerletbodyOutput_S_E_J/Sun.txt", "VelocityVerletbodyOutput_S_E_J/Earth.txt", "VelocityVerletbodyOutput_S_E_J/Jupiter.txt")

eval_planets10("VelocityVerletbodyOutput/Sun.txt", "VelocityVerletbodyOutput/Mercury.txt", "VelocityVerletbodyOutput/Venus.txt", "VelocityVerletbodyOutput/Earth.txt", "VelocityVerletbodyOutput/Mars.txt", "VelocityVerletbodyOutput/Jupiter.txt", "VelocityVerletbodyOutput/Saturn.txt", "VelocityVerletbodyOutput/Uranus.txt", "VelocityVerletbodyOutput/Neptune.txt", "VelocityVerletbodyOutput/Pluto.txt")
