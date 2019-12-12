import sys
import numpy as np
import matplotlib.pyplot as plt

t = []
x = []
y = []
z = []


with open("VVerletoutput.txt") as infile:
    infile.readline()
    lines = infile.readlines()
    for line in lines:
        words = line.split(", ")
        t.append(float(words[0]))
        x.append(float(words[1]))
        y.append(float(words[2]))
        z.append(float(words[3]))

t = np.array(t)
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
