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

dist = np.zeros(len(vec[0]))

times =0

for i in range(len(vec[0])):
    dist[i] = np.sqrt(abs(vec[0][i]**2 + vec[1][i]**2))
    #print(dist[i])
    if(dist[i] < 0.307581): # Crude method, but Mercury's orbital period is 87.969 days. 36600/87.969 is approx. 416 (416.055656), so is hopefully OK
        print(np.arctan(vec[0][i]/vec[1][i]) * 206264.806)
        times += 1

print("Times: " + str(times))
