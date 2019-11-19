import numpy as np
import matplotlib.pyplot as plt
import sys

filename1 = sys.argv[1]
filename2 = sys.argv[2]
filename3 = sys.argv[3]
filename4 = sys.argv[4]

filenames = [filename1, filename2, filename3, filename4]

temperature = [[],[],[],[]]
energy = [[],[],[],[]]
heatcapacity = [[],[],[],[]]
susceptibility = [[],[],[],[]]
meanmagnetization = [[],[],[],[]]

for i in range(len(temperature)):
    with open(filenames[i]) as infile:
        infile.readline()
        lines = infile.readlines()
        for line in lines:
            words = line.split(" , ")
            temperature[i].append(float(words[0]))
            energy[i].append(float(words[1]))
            heatcapacity[i].append(float(words[2]))
            susceptibility[i].append(float(words[3]))
            meanmagnetization[i].append(float(words[4]))

#plt.figure(figsize = (6,10))

for i in range(len(temperature)):
    plt.plot(temperature[i], energy[i], label="$L$ =  %d"%(i*40+20))
plt.title("Strøm som funksjon av spenning for $L = 2 \mu m$")
plt.xlabel("Spenning [V]")
plt.ylabel("Strøm [A]")
plt.legend()
plt.grid()
plt.tight_layout()
#plt.savefig("../bilder/vg-2um-MOSFET.png")
plt.show()

for i in range(len(temperature)):
    plt.plot(temperature[i], heatcapacity[i], label="$L$ =  %d"%(i*40+20))
plt.title("Strøm som funksjon av spenning for $L = 2 \mu m$")
plt.xlabel("Spenning [V]")
plt.ylabel("Strøm [A]")
plt.legend()
plt.grid()
plt.tight_layout()
#plt.savefig("../bilder/vg-2um-MOSFET.png")
plt.show()

for i in range(len(temperature)):
    plt.plot(temperature[i], susceptibility[i], label="$L$ =  %d"%(i*40+20))
plt.title("Strøm som funksjon av spenning for $L = 2 \mu m$")
plt.xlabel("Spenning [V]")
plt.ylabel("Strøm [A]")
plt.legend()
plt.grid()
plt.tight_layout()
#plt.savefig("../bilder/vg-2um-MOSFET.png")
plt.show()

for i in range(len(temperature)):
    plt.plot(temperature[i], meanmagnetization[i], label="$L$ =  %d"%(i*40+20))
plt.title("Strøm som funksjon av spenning for $L = 2 \mu m$")
plt.xlabel("Spenning [V]")
plt.ylabel("Strøm [A]")
plt.legend()
plt.grid()
plt.tight_layout()
#plt.savefig("../bilder/vg-2um-MOSFET.png")
plt.show()
