"""
Last changed: 07.09.2019 23:28 by Erlend
"""

#imports
import numpy as np
import pylab as pl

#create lists for file-data
x = []
v_gen = []
v_elementsspl = []
v_LU = []
u = []

filename = input("What file do you want to read? ") #multiple files useful, we ask for user-input

with open(filename) as infile:
    infile.readline() #discard 1st line
    lines = infile.readlines() #read whole file
    for line in lines:
        words = line.split(", ") #split line per each ", "
        x.append(float(words[0])) #convert and append 1st element
        v_gen.append(float(words[1])) #convert and append 2nd element
        v_spl.append(float(words[2])) #convert and append 3rd element
        v_LU.append(float(words[3])) #convert and append 4th element
        u.append(float(words[4])) #convert and append 5th element

#convert lists to arrays for ease-of-use
x = np.array(x)
v_gen = np.array(v_gen)
v_spl = np.array(v_spl)
v_LU = np.array(v_LU)
u = np.array(u)


"""in hindsight very strange, but it works"""
try:
    n = int(filename[4:-4]) #discard the 4 1st and 4 last characters of filename, (essentially discard "data" and discard ".txt")
except: #if the discardinging in try fails
    print("Could not convert to int. Defaulting to n = 10")
    n = 10 #default value incase of failure


pl.figure(figsize=(12,9)) #change figuresize for better quality image
pl.plot(x,v_gen, label="General Algorithm") #plot of the general algorithm
pl.plot(x,v_spl, label="Special Algorithm") #plot of the special algorithm
pl.plot(x,v_LU, label="LU Algorithm") #plot of the LU-decomp. and solve() method
pl.plot(x,u, label="Exact") #plot "exact" solution as reference
pl.grid()#add grid to plot
pl.title("Thomas-algorithm approximation of Poisson-equation with n = " + str(n)) #add title
pl.xlabel("x") #add x-label
pl.ylabel("u(x)") #add y-label
pl.legend() #enable labels on plot
pl.savefig("data" + str(n) + ".png") #saves plot as image
pl.show() #show plot
