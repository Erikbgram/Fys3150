import sys
import numpy as np
import matplotlib.pyplot as plt


valuelist = []

#input = sys.argv[1]

filename = sys.argv[1]

# Energy
with open("energy/" + filename) as infile:
    infile.readline()
    lines = infile.readlines()
    for line in lines:
        words = line.split(" , ")
        valuelist.append(float(words[0]))

values = np.array(valuelist)
csvalues = np.cumsum(values)
result = np.zeros(len(values))

for i in range(len(values)):
    result[i] = csvalues[i]/(i+1) # cumulative sum / how many cycles have passed

#plt.subplot(2,1,1)
plt.plot(result)
plt.title("Plot of E(MC) for " + filename)
plt.xlabel("MC cycles")
plt.ylabel("Energy")
plt.grid()
plt.savefig("../img/energy_" + filename[:-4] + ".png")
plt.show()
valuelist = []

# Magnet
with open("magnet/" + filename) as infile:
    infile.readline()
    lines = infile.readlines()
    for line in lines:
        words = line.split(" , ")
        valuelist.append(float(words[0]))

values = np.array(valuelist)
csvalues = np.cumsum(values)
result = np.zeros(len(values))

for i in range(len(values)):
    result[i] = csvalues[i]/(i+1) # cumulative sum / how many cycles have passed

#plt.subplot(2,1,2)
plt.plot(result)
plt.title("Plot of M(MC) for " + filename)
plt.xlabel("MC cycles")
plt.ylabel("Magnetization")
plt.grid()
plt.savefig("../img/magnet_" + filename[:-4] + ".png")
plt.show()
