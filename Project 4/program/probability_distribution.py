import sys
import numpy as np
from collections import Counter
import matplotlib.pyplot as plt

complete_energylist = []
energylist = []


energyfile = sys.argv[1]

outputfile = "output" + energyfile[6:]
print(outputfile)

with open(outputfile) as infile: # Reads in the expectation value and std of E
    infile.readline()
    lines = infile.readlines()
    for line in lines:
        words = line.split(" , ")
        avg = float(words[1])
        expect = float(words[2])


with open(sys.argv[1]) as infile: # Reads in the energy values
    infile.readline()
    lines = infile.readlines()
    for line in lines:
        words = line.split()
        energy = float(words[0])
        complete_energylist.append(energy)
        """
        if (energy < (expect+std)) or (energy > (expect-std)): # If E is less than std off from the expectation value
            energylist.append(energy)
        """

complete_energylist = np.array(complete_energylist)
std = complete_energylist.std()

print(avg-std)

for i in range(len(complete_energylist)):
    # skal være AND under
    # as it is now it's broken and helps nothing
    if (complete_energylist[i] < (avg+std)) or (complete_energylist[i] > (avg-std)): # If E is less than std off from the expectation value
        energylist.append(complete_energylist[i])

energy = np.array(energylist)

plt.axhline(y=avg+std, color="r", label="<E>+std")
plt.axhline(y=avg-std, color="g", label="<E>-std")
plt.legend()
plt.plot(complete_energylist)
plt.title("Plot of " + sys.argv[1])
plt.xlabel("iterations")
plt.ylabel("Value")
plt.grid()
plt.show()


plt.hist(energy)
plt.axvline(x=avg+std, color="r", label="<E>+std")
plt.axvline(x=avg-std, color="g", label="<E>-std")
plt.title("Probability of " + sys.argv[1])
plt.xlabel("Energy")
plt.ylabel("P(E)")
plt.show()
