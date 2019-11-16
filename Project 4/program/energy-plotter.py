import matplotlib.pyplot as plt
import numpy as np
import sys

Ei = []



with open(sys.argv[1]) as infile:
    infile.readline()
    lines = infile.readlines()
    for line in lines:
        words = line.split()
        Ei.append(float(words[0]))

Ei = np.array(Ei)

plt.plot(Ei)
plt.title("Energy plot of " + sys.argv[1])
plt.xlabel("iterations")
plt.ylabel("globalE")
plt.grid()
plt.show()
