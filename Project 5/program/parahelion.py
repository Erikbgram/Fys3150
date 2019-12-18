import numpy as np



vec = [[],[]]

with open("../data/VelocityVerlet_S_M_n36600_yr100_newton/Mercury.txt") as infile:
    infile.readline()
    lines = infile.readlines()
    for line in lines:
        words = line.split(" , ")
        vec[0].append(float(words[0]))
        vec[1].append(float(words[1]))

vec = np.array(vec)

for i in range(len(x)):
    print(x[i])
