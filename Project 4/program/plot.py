import sys
import numpy as np
import matplotlib.pyplot as plt


valuelist = []

filename = sys.argv[1]

with open(filename) as infile:
    infile.readline()
    lines = infile.readlines()
    for line in lines:
        words = line.split(" , ")
        valuelist.append(float(words[0]))

values = np.array(valuelist)
csvalues = np.cumsum(values)

result = np.zeros(len(values))

for i in range(1,len(values)):
    result[i] = csvalues[i]/(i) # cumulative sum / how many cycles have passed

print(csvalues)
print(result)


plt.plot(result)
plt.title("Plot of " + sys.argv[1])
plt.xlabel("MC cycles")
plt.ylabel(filename[0:6])
#plt.grid()
plt.ylim(-2.1,-0.8)
plt.show()
