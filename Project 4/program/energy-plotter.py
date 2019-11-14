import matplotlib.pyplot as plt
import numpy as np

Ei = []

with open("energy.txt") as infile:
    infile.readline()
    lines = infile.readlines()
    for line in lines:
        Ei.append(line[0])

Ei = np.array(Ei)

x = np.linspace(0,len(Ei),len(Ei))

plt.plot(x,Ei)
plt.show()
