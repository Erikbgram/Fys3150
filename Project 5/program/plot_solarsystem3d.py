#Sist endret: 17.12.2019 21:46 by Alexandra Jahr Kolstad

import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D     #3D
import matplotlib as mpl    #3D

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
                z[i].append(float(words[2]))

    x = np.array(x)
    y = np.array(y)
    z = np.array(z)


    mpl.rcParams['legend.fontsize'] = 10        #3D, bestemmer størrelsen på legendgreia
    fig = plt.figure(figsize = (10,7))      #3D
    ax = fig.gca(projection='3d')       #3D

    ax.scatter(x[0], y[0], z[0], color = dict["Sun"], label='Sun')
    ax.plot(x[1], y[1], z[1], color = dict["Earth"], label='Earth')
    #ax.plot(x[1], y[1], z[1], color = dict["Mercury"], label='Mercury')
    ax.legend()

    ax.set_facecolor("white")           # Setting the background color
    #plt.title("Two-body solar system with " + filename1[:-22])
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    plt.grid()
    plt.axis("equal")
    plt.tight_layout()

    #title1 = filename1[8] + "_" + filename1[25:-8]
    #title2 = filename1[8] + "_" + filename1[27:-8]
    title1 = filename1[21:24] + "_" + filename1[8] + "_" + filename1[25:-8]
    title2 = filename1[23:26] + "_" + filename1[8] + "_" + filename1[27:-8]

    #print(filename1[28:-14])

    if filename1[8] == "F":
        #plt.title("Two-body solar system with " + filename1[8:20] + " and " + filename1[26:-13] + " iterations.")
        plt.title("Two-body solar system with " + filename1[8:20] + " and " + filename1[25:-8])
        #plt.title("Two-body solar system with " + filename1[8:20] + " and " + filename1[26:-14] + " iterations and 100 years")
        plt.savefig("../img/plot3D_" + title1 + ".png")

    elif filename1[8] == "V":
        #plt.title("Two-body solar system with " + filename1[8:22] + " and " + filename1[28:-13] + " iterations.")
        plt.title("Two-body solar system with " + filename1[8:22] + " and " + filename1[27:-8])
        #plt.title("Two-body solar system with " + filename1[8:22] + " and " + filename1[28:-14] + " iterations and 100 years")
        plt.savefig("../img/plot3D_" + title2 + ".png")

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

    mpl.rcParams['legend.fontsize'] = 10        #3D, bestemme størrelsen på legendgreia
    fig = plt.figure(figsize = (10,7))      #3D
    ax = fig.gca(projection='3d')       #3D

    ax.scatter(x[0], y[0], z[0], color = dict["Sun"], label='Sun')
    ax.plot(x[1], y[1], z[1], color = dict["Earth"], label = "Earth")
    ax.plot(x[2], y[2], z[2], color = dict["Jupiter"],   label = "Jupiter")
    ax.legend()

    ax.set_facecolor("white")           # Setting the background color
    ax.set_xlabel('x')      #få navn på aksene i 3D
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    plt.grid()
    plt.axis("equal")
    plt.tight_layout()

    #title1 = filename1[8] + "_" + filename1[25:-8]
    #title2 = filename1[8] + "_" + filename1[27:-8]
    title1 = filename1[21:26] + "_" + filename1[8] + "_" + filename1[27:-8]
    title2 = filename1[23:28] + "_" + filename1[8] + "_" + filename1[29:-8]

    #print(filename1[21:26])

    if filename1[8] == "F":
        plt.title("Three-body solar system with " + filename1[8:20] + " and " + filename1[27:-8])
        #plt.title("Three-body solar system with " + filename1[8:20] + " and " + filename1[25:-8])
        #plt.title("Three-body solar system with " + filename1[8:20] + " and " + filename1[26:-14] + " iterations and 100 years")
        plt.savefig("../img/plot3D_" + title1 + ".png")

    elif filename1[8] == "V":
        plt.title("Three-body solar system with " + filename1[8:22] + " and " + filename1[29:-8])
        #plt.title("Three-body solar system with " + filename1[8:22] + " and " + filename1[27:-8])
        #plt.title("Three-body solar system with " + filename1[8:22] + " and " + filename1[28:-14] + " iterations and 100 years")
        plt.savefig("../img/plot3D_" + title2 + ".png")


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

    mpl.rcParams['legend.fontsize'] = 10        #3D, bestemme hvor stor størrelse legendgreia skal være
    fig = plt.figure(figsize = (10,7))      #3D
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
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.04),
         ncol=5)
    #ax.legend(loc='center right', bbox_to_anchor=(1, 0.5))


    ax.set_facecolor("white")       # Setting the background color
    plt.title("Ten-body solar system with " + filename1[:-25])
    ax.set_xlabel('x')          #hvordan få labels for x-, y- og z-aksen i 3D
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    plt.grid()
    plt.axis("equal")
    plt.tight_layout()
    plt.savefig("../img/plot3D_10body_" + filename1[:-25] + ".png")
    plt.show()

#--------------------------------------------

#n = "3.00"
#navn = "vel" + n + "pi"

#navn = "b3.0"

m = "10"
n = "8784"
yr = "12"
navn = "m" + m + "_n" + n + "_yr" + yr


#Forward Euler output

#eval_planets2("../data/ForwardEuler_S_E_n" + n + "_yr10/Sun.txt", "../data/ForwardEuler_S_E_n" + n + "_yr10/Earth.txt")

#eval_planets2("../data/ForwardEuler_S_E_" + navn + "/Sun.txt", "../data/ForwardEuler_S_E_" + navn + "/Earth.txt")

#eval_planets2("../data/ForwardEuler_S_M_" + navn + "/Sun.txt", "../data/ForwardEuler_S_M_" + navn + "/Mercury.txt")

eval_planets3("../data/ForwardEuler_S_E_J_" + navn + "/Sun.txt", "../data/ForwardEuler_S_E_J_" + navn + "/Earth.txt", "../data/ForwardEuler_S_E_J_" + navn + "/Jupiter.txt")

#eval_planets10("../data/ForwardEuler_10body_n3660_yr10/Sun.txt", "../data/ForwardEuler_10body_n3660_yr10/Mercury.txt", "../data/ForwardEuler_10body_n3660_yr10/Venus.txt", "../data/ForwardEuler_10body_n3660_yr10/Earth.txt", "../data/ForwardEuler_10body_n3660_yr10/Mars.txt", "../data/ForwardEuler_10body_n3660_yr10/Jupiter.txt", "../data/ForwardEuler_10body_n3660_yr10/Saturn.txt", "../data/ForwardEuler_10body_n3660_yr10/Uranus.txt", "../data/ForwardEuler_10body_n3660_yr10/Neptune.txt", "../data/ForwardEuler_10body_n3660_yr10/Pluto.txt")





#Velocity Verlet output

#eval_planets2("../data/VelocityVerlet_S_E_n" + n + "_yr10/Sun.txt", "../data/VelocityVerlet_S_E_n" + n + "_yr10/Earth.txt")

#eval_planets2("../data/VelocityVerlet_S_E_" + navn + "/Sun.txt", "../data/VelocityVerlet_S_E_" + navn + "/Earth.txt")

#eval_planets2("../data/VelocityVerlet_S_M_" + navn + "/Sun.txt", "../data/VelocityVerlet_S_M_" + navn + "/Mercury.txt")

eval_planets3("../data/VelocityVerlet_S_E_J_" + navn + "/Sun.txt", "../data/VelocityVerlet_S_E_J_" + navn + "/Earth.txt", "../data/VelocityVerlet_S_E_J_" + navn + "/Jupiter.txt")

#eval_planets10("../data/VelocityVerlet/Sun.txt", "../data/VelocityVerlet/Mercury.txt", "../data/VelocityVerlet/Venus.txt", "../data/VelocityVerlet/Earth.txt", "../data/VelocityVerlet/Mars.txt", "../data/VelocityVerlet/Jupiter.txt", "../data/VelocityVerlet/Saturn.txt", "../data/VelocityVerlet/Uranus.txt", "../data/VelocityVerlet/Neptune.txt", "../data/VelocityVerlet/Pluto.txt")
