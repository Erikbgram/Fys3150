import sys
import numpy as np
import matplotlib.pyplot as plt


valuelist = []

filename = sys.argv[1] # Must input accept/NAME OF FILE


with open("accept/" + filename) as infile:
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


plt.plot(result)
plt.title("Plot of Accepted spins as a function of MC-cycles for " + filename)
plt.xlabel("MC cycles")
plt.ylabel("Accepted Spins")
plt.grid()
plt.savefig("../img/accept_" + filename[:-4] + ".png")
plt.show()
