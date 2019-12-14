import sys
import numpy as np
import matplotlib.pyplot as plt

x = []
y = []
z = []


with open("BodyOutput/Earth.txt") as infile:
    infile.readline()
    lines = infile.readlines()
    for line in lines:
        words = line.split(", ")
        x.append(float(words[0]))
        y.append(float(words[1]))
        z.append(float(words[2]))

x = np.array(x)
y = np.array(y)
z = np.array(z)

plt.plot(x,y)
plt.plot(0, 0, "yo")
plt.title("Earth")
plt.xlabel("x")
plt.ylabel("y")
plt.grid()
#plt.tight_layout()
#plt.savefig("../img/earth.png")
plt.show()
