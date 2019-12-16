#Sist endret: 16.12.2019 14:19 by Alexandra Jahr Kolstad

import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D     #3D
import matplotlib as mpl    #3D

dict = {"Sun":"#ffdf22" , "Mercury":"#d3d3d3" , "Venus":"#cde7ff" , "Earth":"#197619" , "Mars":"#cf6e28" , "Jupiter":"#e69500" , "Saturn":"#b9a17e" , "Uranus":"#78d6fb" , "Neptune":"#1414fb" , "Pluto":"#ddb28d"}

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
                z[i].append(float(words[2]))

    x = np.array(x)
    y = np.array(y)
    z = np.array(z)


    mpl.rcParams['legend.fontsize'] = 10        #3D, bestemmer størrelsen på legendgreia
    fig = plt.figure()      #3D
    ax = fig.gca(projection='3d')       #3D

    ax.scatter(x[0], y[0], z[0], color = dict["Sun"], label='Sun')
    ax.plot(x[1], y[1], z[1], color = dict["Earth"], label='Earth')
    ax.legend()

    ax.set_facecolor("white")           # Setting the background color
    plt.title("Two-body solar system with " + filename1[:-22])
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    plt.grid()
    plt.axis("equal")
    plt.tight_layout()
    #plt.savefig("../img/earth.png")
    plt.show()

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

    mpl.rcParams['legend.fontsize'] = 7        #3D, bestemme størrelsen på legendgreia
    fig = plt.figure()      #3D
    ax = fig.gca(projection='3d')       #3D

    ax.scatter(x[0], y[0], z[0], color = dict["Sun"], label='Sun')
    ax.plot(x[1], y[1], z[1], color = dict["Earth"], label = "Earth")
    ax.plot(x[2], y[2], z[2], color = dict["Jupiter"],   label = "Jupiter")
    ax.legend()

    ax.set_facecolor("white")           # Setting the background color
    plt.title("Three-body solar system with " + filename1[:-24])
    ax.set_xlabel('x')      #få navn på aksene i 3D
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    plt.grid()
    plt.axis("equal")
    plt.tight_layout()
    #plt.savefig("../img/earth.png")
    plt.show()

def eval_planets10(filename1, filename2, filename3, filename4, filename5, filename6, filename7, filename8, filename9, filename10):

    filenames = [filename1, filename2, filename3, filename4, filename5, filename6, filename7, filename8, filename9, filename10]

    x = [[],[],[],[],[],[],[],[],[],[]]
    y = [[],[],[],[],[],[],[],[],[],[]]
    z = [[],[],[],[],[],[],[],[],[],[]]

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

    mpl.rcParams['legend.fontsize'] = 7        #3D, bestemme hvor stor størrelse legendgreia skal være
    fig = plt.figure()      #3D
    ax = fig.gca(projection='3d')       #3D

    ax.scatter(x[0], y[0], z[0], color = dict["Sun"], label='Sun')
    ax.plot(x[1], y[1], z[1], color = dict["Mercury"], label = "Mercury")
    ax.plot(x[2], y[2], z[2], color = dict["Venus"],   label = "Venus")
    ax.plot(x[3], y[3], z[3], color = dict["Earth"],   label = "Earth")
    ax.plot(x[4], y[4], z[4], color = dict["Mars"],    label = "Mars")
    ax.plot(x[5], y[5], z[5], color = dict["Jupiter"], label = "Jupiter")
    ax.plot(x[6], y[6], z[6], color = dict["Saturn"],  label = "Saturn")
    ax.plot(x[7], y[7], z[7], color = dict["Uranus"],  label = "Uranus")
    ax.plot(x[8], y[8], z[8], color = dict["Neptune"], label = "Neptune")
    ax.plot(x[9], y[9], z[9], color = dict["Pluto"],   label = "Pluto")
    ax.legend(loc = "best")

    ax.set_facecolor("white")       # Setting the background color
    plt.title("Ten-body solar system with " + filename1[:-18])
    ax.set_xlabel('x')          #hvordan få labels for x-, y- og z-aksen i 3D
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    plt.grid()
    plt.axis("equal")
    plt.tight_layout()
    #plt.savefig("../img/earth.png")
    plt.show()

def main():

    #Forward Euler output

    eval_planets2("ForwardEulerbodyOutput_S_E/Sun.txt", "ForwardEulerbodyOutput_S_E/Earth.txt")

    #eval_planets3("ForwardEulerbodyOutput_S_E_J/Sun.txt", "ForwardEulerbodyOutput_S_E_J/Earth.txt", "ForwardEulerbodyOutput_S_E_J/Jupiter.txt")

    eval_planets10("ForwardEulerbodyOutput/Sun.txt", "ForwardEulerbodyOutput/Mercury.txt", "ForwardEulerbodyOutput/Venus.txt", "ForwardEulerbodyOutput/Earth.txt", "ForwardEulerbodyOutput/Mars.txt", "ForwardEulerbodyOutput/Jupiter.txt", "ForwardEulerbodyOutput/Saturn.txt", "ForwardEulerbodyOutput/Uranus.txt", "ForwardEulerbodyOutput/Neptune.txt", "ForwardEulerbodyOutput/Pluto.txt")

    #Velocity verlet output

    #eval_planets2("VelocityVerletbodyOutput_S_E/Sun.txt", "VelocityVerletbodyOutput_S_E/Earth.txt")

    #eval_planets3(""VelocityVerletbodyOutput_S_E_J/Sun.txt", "VelocityVerletbodyOutput_S_E_J/Earth.txt", "VelocityVerletbodyOutput_S_E_J/Jupiter.txt")

    eval_planets10("VelocityVerletbodyOutput/Sun.txt", "VelocityVerletbodyOutput/Mercury.txt", "VelocityVerletbodyOutput/Venus.txt", "VelocityVerletbodyOutput/Earth.txt", "VelocityVerletbodyOutput/Mars.txt", "VelocityVerletbodyOutput/Jupiter.txt", "VelocityVerletbodyOutput/Saturn.txt", "VelocityVerletbodyOutput/Uranus.txt", "VelocityVerletbodyOutput/Neptune.txt", "VelocityVerletbodyOutput/Pluto.txt")

if __name__ == "__main__":
    main()
