import sys
import numpy as np
from collections import Counter
import matplotlib.pyplot as plt

file = sys.argv[1]

energy = []
magnet = []

# Energy
with open("energy/" + file) as infile: # Reads in the energy values
    infile.readline()
    lines = infile.readlines()
    for line in lines:
        words = line.split()
        energy.append(float(words[0]))

# Magnetization
with open("magnet/" + file) as infile: # Reads in the energy values
    infile.readline()
    lines = infile.readlines()
    for line in lines:
        words = line.split()
        magnet.append(float(words[0]))



energy = np.array(energy)
estd = energy.std()
emean = energy.mean()
evar = np.var(energy)

magnet = np.array(magnet)
mstd = magnet.std()
mmean = magnet.mean()
mvar = np.var(magnet)


"""
plt.axhline(y=emean-estd, color="r")
plt.axhline(y=emean, color="g", label="mean")
plt.axhline(y=emean+estd, color="r", label="std")
plt.plot(energy)
plt.title("Energy plot of " + sys.argv[1])
plt.xlabel("Energy")
plt.ylabel("P(E)")
plt.grid()
plt.show()
"""
# Energy
plt.axvline(emean-estd, linestyle="--", color="r")
plt.axvline(emean, linestyle="--", color="g", label="mean")
plt.axvline(emean+estd, linestyle="--", color="r", label="std")
plt.hist(energy[10000:], bins=114, density=True, stacked=True) # Tried various bins, but this gave the fewest artifacts
plt.title("Energy Probability distribution of \n" + sys.argv[1])
plt.xlabel("Energy")
plt.ylabel("P(E)")
plt.legend()
plt.tight_layout()
plt.savefig("../img/energyhistogram_" + file[:-11] + file[-10:-4] + ".png")
#plt.savefig("../img/energyhistogram_" + sys.argv[1][:-4] + ".png")
plt.show()


# Magnetization
plt.axvline(mmean-mstd, linestyle="--", color="r")
plt.axvline(mmean, linestyle="--", color="g", label="mean")
plt.axvline(mmean+mstd, linestyle="--", color="r", label="std")
plt.hist(magnet[10000:], bins=126, density=True, stacked=True) # Tried various bins, but this gave the fewest artifacts
plt.title("Magnetization Probability distribution of \n" + sys.argv[1])
plt.xlabel("Magnetization")
plt.ylabel("P(E)")
plt.legend()
plt.tight_layout()
plt.savefig("../img/magnethistogram_" + file[:-11] + file[-10:-4] + ".png")
#plt.savefig("../img/magnethistogram_" + sys.argv[1][:-4] + ".png")
plt.show()
