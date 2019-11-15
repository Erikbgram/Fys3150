import matplotlib.pyplot as plt
import numpy as np

Ei = []

with open("energy.txt") as infile:
    infile.readline()
    lines = infile.readlines()
    for line in lines:
        words = line.split()
        Ei.append(float(words[0]))

Ei = np.array(Ei)

plt.plot(Ei)
plt.xlabel("iterations")
plt.ylabel("globalE")
plt.grid()
plt.show()
