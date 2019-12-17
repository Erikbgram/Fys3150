#Sist endret: 17.12.2019 21:46 by Alexandra Jahr Kolstad

import sys
import numpy as np
import matplotlib.pyplot as plt

dict = {"Sun":"#ffdf22" , "Mercury":"#d3d3d3" , "Venus":"#c04e01" , "Earth":"#197619" , "Mars":"#cf6e28" , "Jupiter":"#e69500" , "Saturn":"#b9a17e" , "Uranus":"#78d6fb" , "Neptune":"#1414fb" , "Pluto":"#ddb28d"}

def eval_planets2(filename1, filename2):

    filenames = [filename1, filename2]

    x = [[],[]]
    y = [[],[]]

    for i in range(len(x)):
        with open(filenames[i]) as infile:
            infile.readline()
            lines = infile.readlines()
            for line in lines:
                words = line.split(" , ")
                x[i].append(float(words[0]))
                y[i].append(float(words[1]))

    x = np.array(x)
    y = np.array(y)

    plt.figure(figsize=(10,6))
    ax = plt.axes()
    ax.set_facecolor("grey")           # Setting the background color
    plt.scatter(x[0], y[0], color = dict["Sun"], label = "Sun")
    plt.plot(x[1],y[1], color = dict["Earth"], label = "Earth")

    plt.legend()
    plt.title("Two-body solar system with " + filename1[:-22])
    plt.xlabel("x")
    plt.ylabel("y")
    plt.grid()
    plt.axis("equal")
    plt.tight_layout()
    #plt.savefig("../img/earth.png")
    plt.show()

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

    x = np.array(x)
    y = np.array(y)

    plt.figure(figsize=(10,6))
    ax = plt.axes()
    ax.set_facecolor("grey")            # Setting the background color
    plt.scatter(x[0], y[0], color = dict["Sun"], label = "Sun")
    plt.plot(x[1],y[1], color = dict["Earth"], label = "Earth")
    plt.plot(x[2],y[2], color = dict["Jupiter"], label = "Jupiter")

    plt.legend()
    plt.title("Three-body solar system with " + filename1[:-24])
    plt.xlabel("x")
    plt.ylabel("y")
    plt.grid()
    plt.axis("equal")
    plt.tight_layout()
    #plt.savefig("../img/earth.png")
    plt.show()

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

    x = np.array(x)
    y = np.array(y)

    plt.figure(figsize=(10,6))
    ax = plt.axes()
    ax.set_facecolor("grey")    # Setting the background color

    plt.scatter(x[0], y[0], color = dict["Sun"], label = "Sun")
    plt.plot(x[1],y[1], color = dict["Mercury"], label = "Mercury")
    plt.plot(x[2],y[2], color = dict["Venus"], label = "Venus")
    plt.plot(x[3],y[3], color = dict["Earth"], label = "Earth")
    plt.plot(x[4],y[4], color = dict["Mars"], label = "Mars")
    plt.plot(x[5],y[5], color = dict["Jupiter"], label = "Jupiter")
    plt.plot(x[6],y[6], color = dict["Saturn"], label = "Saturn")
    plt.plot(x[7],y[7], color = dict["Uranus"], label = "Uranus")
    plt.plot(x[8],y[8], color = dict["Neptune"], label = "Neptune")
    plt.plot(x[9],y[9], color = dict["Pluto"], label = "Pluto")
    plt.legend()
    plt.title("Ten-body solar system with " + filename1[:-18])
    plt.xlabel("x")
    plt.ylabel("y")
    plt.grid()
    plt.axis("equal")
    plt.tight_layout()
    #plt.savefig("../img/earth.png")
    plt.show()

#----------------------------------------


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
