import matplotlib.pyplot as plt
import numpy as np
import sys

values = []

file = sys.argv[1]

if file[0:6] == "energy":

    outputfile = "output/" + file[7:]

    with open(outputfile) as infile: # Reads in the expectation value and std of E
        infile.readline()
        lines = infile.readlines()
        for line in lines:
            words = line.split(" , ")
            avg = float(words[1])
            expect = float(words[2])

    with open(sys.argv[1]) as infile:
        infile.readline()
        lines = infile.readlines()
        for line in lines:
            words = line.split()
            values.append(float(words[0]))

    values = np.array(values)

    std = values.std()

    plt.axhline(y=avg+std, color="r", label="<E>+std")
    plt.axhline(y=avg-std, color="g", label="<E>-std")
    plt.legend()

else:
    with open(sys.argv[1]) as infile:
        infile.readline()
        lines = infile.readlines()
        for line in lines:
            words = line.split()
            values.append(float(words[0]))

    values = np.array(values)

plt.plot(np.cumsum(values))
plt.xscale("log")
plt.title("Plot of " + sys.argv[1])
plt.xlabel("MC cycles")
plt.ylabel("Value")
plt.grid()
plt.show()
