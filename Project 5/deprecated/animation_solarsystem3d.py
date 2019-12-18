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

    #mpl.rcParams['animation.convert_path'] = r'C:\Program Files\ImageMagick\convert'
    #mpl.rcParams['animation.ffmpeg_path'] = r'C:\Program Files\ffmpeg\bin\ffmpeg.exe'


    line1, = ax.plot(x[0], y[0], z[0], color = dict["Sun"], label = "Sun")
    line2, = ax.plot(x[1], y[1], z[1], color = dict["Earth"], label = "Earth")

    #ani = animation.FuncAnimation(fig, update_planets2, frames = len(x[0]), fargs=[x, y, z, line1, line2],
    #                  interval=10, blit=True)

    # Set up formatting for the movie files
    #Writer = animation.writers['ffmpeg']
    #writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
    #writer = Writer(fps=15, bitrate=1800)

    ax.legend()
    ax.set_facecolor("white")           # Setting the background color
    #plt.title("Two-body solar system with " + filename1[:-22])
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    plt.grid()
    plt.axis("equal")
    plt.tight_layout()

    #ani = animation.FuncAnimation(fig, update_planets2, frames = len(x[0]), fargs=[x, y, z, line1, line2],
                      #interval=1, blit=True)

    ani = animation.FuncAnimation(fig, update_planets2, frames = len(x[0]), fargs=[x, y, z, line1, line2],
                      interval=1, blit=True)

    #mywriter = animation.FFMpegWriter()

    #ani.save("ani3D_2body_" + filename1[:-22] + ".png", writer=writer)
    #ani.save("ani3D_2body_" + filename1[:-22] + ".gif", writer = "imagemagick")

    if filename1[0] == "F":
        plt.title("Two-body solar system with " + filename1[:12] + " and " + filename1[18:-10] + " iterations.")
        ani.save("../img/ani3D_2body_" + filename1[:12] + "_" + filename1[18:-10] + ".gif", writer = "imagemagick", fps = 60)
    elif filename1[0] == "V":
        plt.title("Two-body solar system with " + filename1[:13] + " and " + filename1[18:5] + " iterations.")
        #ani.save("ani3D_2body.gif", writer = "imagemagick")

    plt.show()

def update_planets2(num, x, y, z, line1, line2):
    line1.set_data(x[0][:num], y[0][:num])
    line2.set_data(x[1][:num], y[1][:num])
    line1.set_3d_properties(z[0][:num])
    line2.set_3d_properties(z[1][:num])
    return [line1,line2]

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

    ani = animation.FuncAnimation(fig, update_planets10, frames = len(x[0]), fargs=[x, y, z, line1, line2, line3, line4, line5, line6, line7, line8, line9, line10],
                      interval=1, blit=True)

    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.04),
         ncol=5)
    ax.set_facecolor("white")       # Setting the background color
    plt.title("Ten-body solar system with " + filename1[:-25])
    ax.set_xlabel('x')          #hvordan få labels for x-, y- og z-aksen i 3D
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    plt.grid()
    plt.axis("equal")
    plt.tight_layout()
    #plt.savefig("../img/plot3D_10body_" + filename1[:-25] + ".png")
    plt.show()

def update_planets10(num, x, y, z, line1, line2, line3, line4, line5, line6, line7, line8, line9, line10):
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
    line1.set_3d_properties(z[0][:num])
    line2.set_3d_properties(z[1][:num])
    line3.set_3d_properties(z[2][:num])
    line4.set_3d_properties(z[3][:num])
    line5.set_3d_properties(z[4][:num])
    line6.set_3d_properties(z[5][:num])
    line7.set_3d_properties(z[6][:num])
    line8.set_3d_properties(z[7][:num])
    line9.set_3d_properties(z[8][:num])
    line10.set_3d_properties(z[9][:num])
    return [line1,line2,line3,line4,line5,line6,line7,line8,line9,line10]


#-----------------------------------------------------------------------------------------------

#Forward Euler output

eval_planets2("ForwardEulerbodyOutput_S_E/Sun.txt", "ForwardEulerbodyOutput_S_E/Earth.txt")

#eval_planets3("ForwardEulerbodyOutput_S_E_J/Sun.txt", "ForwardEulerbodyOutput_S_E_J/Earth.txt", "ForwardEulerbodyOutput_S_E_J/Jupiter.txt")

#eval_planets10("ForwardEulerbodyOutput_10body/Sun.txt", "ForwardEulerbodyOutput_10body/Mercury.txt", "ForwardEulerbodyOutput_10body/Venus.txt", "ForwardEulerbodyOutput_10body/Earth.txt", "ForwardEulerbodyOutput_10body/Mars.txt", "ForwardEulerbodyOutput_10body/Jupiter.txt", "ForwardEulerbodyOutput_10body/Saturn.txt", "ForwardEulerbodyOutput_10body/Uranus.txt", "ForwardEulerbodyOutput_10body/Neptune.txt", "ForwardEulerbodyOutput_10body/Pluto.txt")

#Velocity verlet output

#eval_planets2("VelocityVerletbodyOutput_S_E/Sun.txt", "VelocityVerletbodyOutput_S_E/Earth.txt")

#eval_planets3(""VelocityVerletbodyOutput_S_E_J/Sun.txt", "VelocityVerletbodyOutput_S_E_J/Earth.txt", "VelocityVerletbodyOutput_S_E_J/Jupiter.txt")

#eval_planets10("VelocityVerletbodyOutput_10body/Sun.txt", "VelocityVerletbodyOutput_10body/Mercury.txt", "VelocityVerletbodyOutput_10body/Venus.txt", "VelocityVerletbodyOutput_10body/Earth.txt", "VelocityVerletbodyOutput_10body/Mars.txt", "VelocityVerletbodyOutput_10body/Jupiter.txt", "VelocityVerletbodyOutput_10body/Saturn.txt", "VelocityVerletbodyOutput_10body/Uranus.txt", "VelocityVerletbodyOutput_10body/Neptune.txt", "VelocityVerletbodyOutput_10body/Pluto.txt")
