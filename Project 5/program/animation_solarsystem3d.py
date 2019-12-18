#Sist endret: 17.12.2019 21:45 by Alexandra Jahr Kolstad

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D     #3D
import matplotlib as mpl    #3D
import matplotlib.animation as animation
#mpl.use("Agg")

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


    mpl.rcParams['legend.fontsize'] = 10        #3D, bestemmer størrelsen på legendgreia
    fig = plt.figure(figsize = (10,7))      #3D
    ax = fig.gca(projection='3d')       #3D

    ax.scatter(x[0],y[0],z[0], color = dict["Sun"], label = "Sun" )

    #line1, = ax.plot(x[0], y[0], z[0], color = dict["Sun"], label = "Sun")
    line2, = ax.plot(x[1], y[1], z[1], color = dict["Earth"], label = "Earth")

    ax.set_facecolor("white")           # Setting the background color
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    plt.grid()
    plt.axis("equal")
    plt.tight_layout()

    ax.legend()

    #ani = animation.FuncAnimation(fig, update_planets2, frames = len(x[0]), fargs=[x, y, z, line1, line2],
                      #interval=1, blit=True)

    ani = animation.FuncAnimation(fig, update_planets2, frames = len(x[0]), fargs=[x, y, z, line2],
                      interval=1, blit=True)



    title1 = filename1[8] + "_520"

    #Forward Euler
    #plt.title("Two-body solar system with " + filename1[8:20] + " and " + filename1[26:-13] + " iterations.")
    #ani.save("../img/ani3D_2body_" + title1 + ".gif", writer = "imagemagick", fps = 60)

    #Velocity Verlet
    #plt.title("Two-body solar system with " + filename1[8:22] + " and " + filename1[28:-13] + " iterations.")
    #ani.save("../img/ani3D_2body_" + title1 + ".gif", writer = "imagemagick", fps = 60)

    if filename1[8] == "F":
        #plt.title("Two-body solar system with " + filename1[:12] + " and " + filename1[18:-10] + " iterations.")
        #ani.save("../img/ani3D_2body_" + filename1[:12] + "_" + filename1[18:-10] + ".gif", writer = "imagemagick", fps = 60)

        plt.title("Two-body solar system with " + filename1[8:20] + " and " + filename1[26:-13] + " iterations.")
        ani.save("../img/ani3D_2body_" + title1 + ".gif", writer = "imagemagick", fps = 60)

    elif filename1[8] == "V":
        #plt.title("Two-body solar system with " + filename1[:13] + " and " + filename1[18:5] + " iterations.")
        #ani.save("ani3D_2body.gif", writer = "imagemagick")

        plt.title("Two-body solar system with " + filename1[8:22] + " and " + filename1[28:-13] + " iterations.")
        ani.save("../img/ani3D_2body_" + title1 + ".gif", writer = "imagemagick", fps = 60)


    #plt.show()

def update_planets2(num, x, y, z, line2):
    #line1.set_data(x[0][:num*2], y[0][:num*2])
    #line2.set_data(x[1][:num*2], y[1][:num*2])
    #line1.set_3d_properties(z[0][:num*2])
    #line2.set_3d_properties(z[1][:num*2])

    #line1.set_data(x[0][:num], y[0][:num])
    line2.set_data(x[1][:num], y[1][:num])
    #line1.set_3d_properties(z[0][:num])
    line2.set_3d_properties(z[1][:num])
    #return [line1,line2]
    return [line2]

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


    mpl.rcParams['legend.fontsize'] = 10        #3D, bestemmer størrelsen på legendgreia
    fig = plt.figure(figsize = (10,7))      #3D
    ax = fig.gca(projection='3d')       #3D

    line1, = plt.plot(x[0], y[0], z[0], color = dict["Sun"], label = "Sun")
    line2, = plt.plot(x[1], y[1], z[1], color = dict["Earth"], label = "Earth")
    line3, = plt.plot(x[2], y[2], z[2], color = dict["Jupiter"], label = "Jupiter")

    ani = animation.FuncAnimation(fig, update_planets3, frames = len(x[0]), fargs=[x, y, z, line1, line2, line3],
                      interval=1, blit=True)

    ax.legend()
    ax.set_facecolor("white")           # Setting the background color
    plt.title("Three-body solar system with " + filename1[:-24])
    ax.set_xlabel('x')      #få navn på aksene i 3D
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    plt.grid()
    plt.axis("equal")
    plt.tight_layout()
    #plt.savefig("../img/plot3D_3body_" + filename1[:-24] + ".png")
    plt.show()

def update_planets3(num, x, y, z, line1, line2, line3):
    line1.set_data(x[0][:num], y[0][:num])
    line2.set_data(x[1][:num], y[1][:num])
    line3.set_data(x[2][:num], y[2][:num])
    line1.set_3d_properties(z[0][:num])
    line2.set_3d_properties(z[1][:num])
    line3.set_3d_properties(z[2][:num])
    return [line1,line2,line3]

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

    mpl.rcParams['legend.fontsize'] = 10        #3D, bestemmer størrelsen på legendgreia
    fig = plt.figure(figsize = (10,7))      #3D
    ax = fig.gca(projection='3d')       #3D

    line1, = plt.plot(x[0], y[0], z[0], color = dict["Sun"], label = "Sun")
    line2, = plt.plot(x[1], y[1], z[1],  color = dict["Mercury"], label = "Mercury")
    line3, = plt.plot(x[2], y[2], z[2], color = dict["Venus"], label = "Venus")
    line4, = plt.plot(x[3], y[3], z[3], color = dict["Earth"], label = "Earth")
    line5, = plt.plot(x[4], y[4], z[4], color = dict["Mars"], label = "Mars")
    line6, = plt.plot(x[5], y[5], z[5], color = dict["Jupiter"], label = "Jupiter")
    line7, = plt.plot(x[6], y[6], z[6], color = dict["Saturn"], label = "Saturn")
    line8, = plt.plot(x[7], y[7], z[7], color = dict["Uranus"], label = "Uranus")
    line9, = plt.plot(x[8], y[8], z[8], color = dict["Neptune"], label = "Neptune")
    line10, = plt.plot(x[9], y[9], z[9], color = dict["Pluto"], label = "Pluto")


    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.04),
         ncol=5)
    ax.set_facecolor("white")       # Setting the background color
    #plt.title("Ten-body solar system with " + filename1[:-25])
    ax.set_xlabel('x')          #hvordan få labels for x-, y- og z-aksen i 3D
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    plt.grid()
    plt.axis("equal")
    plt.tight_layout()

    ani = animation.FuncAnimation(fig, update_planets10, frames = int(len(x[0])/2), fargs=[x, y, z, line1, line2, line3, line4, line5, line6, line7, line8, line9, line10],
                      interval=1, blit=True)

    plt.title("Ten-body solar system with " + filename1[8:-8] + " and 3660 iterations.")
    ani.save("../img/ani3D_10body_3660" + filename1[8:-8] + ".mp4", writer = "imagemagick", fps = 60)

    plt.show()

def update_planets10(num, x, y, z, line1, line2, line3, line4, line5, line6, line7, line8, line9, line10):
    line1.set_data(x[0][:num*2], y[0][:num*2])
    line2.set_data(x[1][:num*2], y[1][:num*2])
    line3.set_data(x[2][:num*2], y[2][:num*2])
    line4.set_data(x[3][:num*2], y[3][:num*2])
    line5.set_data(x[4][:num*2], y[4][:num*2])
    line6.set_data(x[5][:num*2], y[5][:num*2])
    line7.set_data(x[6][:num*2], y[6][:num*2])
    line8.set_data(x[7][:num*2], y[7][:num*2])
    line9.set_data(x[8][:num*2], y[8][:num*2])
    line10.set_data(x[9][:num*2], y[9][:num*2])
    line1.set_3d_properties(z[0][:num*2])
    line2.set_3d_properties(z[1][:num*2])
    line3.set_3d_properties(z[2][:num*2])
    line4.set_3d_properties(z[3][:num*2])
    line5.set_3d_properties(z[4][:num*2])
    line6.set_3d_properties(z[5][:num*2])
    line7.set_3d_properties(z[6][:num*2])
    line8.set_3d_properties(z[7][:num*2])
    line9.set_3d_properties(z[8][:num*2])
    line10.set_3d_properties(z[9][:num*2])
    return [line1,line2,line3,line4,line5,line6,line7,line8,line9,line10]


#-----------------------------------------------------------------------------------------------

n = "520"

#Forward Euler output

eval_planets2("../data/ForwardEuler_S_E_n" + n + "_yr10/Sun.txt", "../data/ForwardEuler_S_E_n" + n + "_yr10/Earth.txt")

#eval_planets3("ForwardEulerbodyOutput_S_E_J/Sun.txt", "ForwardEulerbodyOutput_S_E_J/Earth.txt", "ForwardEulerbodyOutput_S_E_J/Jupiter.txt")

#eval_planets10("../data/ForwardEuler_10body_n3660_yr10/Sun.txt", "../data/ForwardEuler_10body_n3660_yr10/Mercury.txt", "../data/ForwardEuler_10body_n3660_yr10/Venus.txt", "../data/ForwardEuler_10body_n3660_yr10/Earth.txt", "../data/ForwardEuler_10body_n3660_yr10/Mars.txt", "../data/ForwardEuler_10body_n3660_yr10/Jupiter.txt", "../data/ForwardEuler_10body_n3660_yr10/Saturn.txt", "../data/ForwardEuler_10body_n3660_yr10/Uranus.txt", "../data/ForwardEuler_10body_n3660_yr10/Neptune.txt", "../data/ForwardEuler_10body_n3660_yr10/Pluto.txt")

#Velocity Verlet output

eval_planets2("../data/VelocityVerlet_S_E_n" + n + "_yr10/Sun.txt", "../data/VelocityVerlet_S_E_n" + n + "_yr10/Earth.txt")

#eval_planets3(""VelocityVerletbodyOutput_S_E_J/Sun.txt", "VelocityVerletbodyOutput_S_E_J/Earth.txt", "VelocityVerletbodyOutput_S_E_J/Jupiter.txt")

#eval_planets10("../data/VelocityVerlet/Sun.txt", "../data/VelocityVerlet/Mercury.txt", "../data/VelocityVerlet/Venus.txt", "../data/VelocityVerlet/Earth.txt", "../data/VelocityVerlet/Mars.txt", "../data/VelocityVerlet/Jupiter.txt", "../data/VelocityVerlet/Saturn.txt", "../data/VelocityVerlet/Uranus.txt", "../data/VelocityVerlet/Neptune.txt", "../data/VelocityVerlet/Pluto.txt")
