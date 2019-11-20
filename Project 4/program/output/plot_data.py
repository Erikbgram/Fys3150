# Sist endret: 20.11.2019 14:17 av Alexandra Jahr Kolstad

import numpy as np
import matplotlib.pyplot as plt
import sys


filename1 = sys.argv[1] + "_n1000000_Tvar_2.txt"
filename2 = sys.argv[2] + "_n1000000_Tvar_2.txt"
filename3 = sys.argv[3] + "_n1000000_Tvar_2.txt"
filename4 = sys.argv[4] + "_n1000000_Tvar_2.txt"


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

#plt.figure(figsize = (8,10))

for i in range(len(temperature)):
    plt.plot(temperature[i], energy[i], label = "$L$ = %d" %(40+i*20))
plt.title("Energy as a function of temperature.")
plt.xlabel("T [energy]")
plt.ylabel("E [dimless]")
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig("../../img/tempvsenergy2.png")
plt.show()

for i in range(len(temperature)):
    plt.plot(temperature[i], heatcapacity[i], label = "$L$ = %d" %(40+i*20))
plt.title("Specific heat as a function of temperature.")
plt.xlabel("T [energy]")
plt.ylabel("C$_V$ [energy/temperature]")
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig("../../img/tempvsspecificheat2.png")
plt.show()

for i in range(len(temperature)):
    plt.plot(temperature[i], susceptibility[i], label = "$L$ = %d" %(40+i*20))
plt.title("Susceptibility as a function of temperature.")
plt.xlabel("T [energy]")
plt.ylabel("$\chi$ [1/energy]")
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig("../../img/tempvssusceptibility2.png")
plt.show()

for i in range(len(temperature)):
    plt.plot(temperature[i], meanmagnetization[i], label = "$L$ = %d" %(40+i*20))
plt.title("Mean magnetization as a function of temperature.")
plt.xlabel("T [energy]")
plt.ylabel("|M| [dimless]")
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig("../../img/tempvsmeanmagnetization2.png")
plt.show()
