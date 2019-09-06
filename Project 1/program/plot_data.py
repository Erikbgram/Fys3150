"""
Last changed: 06.09.2019 16:32 by Erlend
"""

import numpy as np
import pylab as pl

x = []
v = []
u = []

filename = input("What file do you want to read? ")

with open(filename) as infile:
    infile.readline()
    lines = infile.readlines()
    for line in lines:
        words = line.split(", ")
        x.append(float(words[0]))
        v.append(float(words[1]))
        u.append(float(words[2]))

x = np.array(x)
v = np.array(v)
u = np.array(u)

try:
    n = int(filename[12:-4])
except:
    print("Could not convert to int. Defaulting to n = 10")
    n = 10


pl.plot(x,v, label="Approximation")
pl.plot(x,u, label="Exact")
pl.grid()
pl.title("Thomas-Algorithm approximation of poisson-equation with n = " + str(n))
pl.xlabel("x")
pl.ylabel("u(x)")
pl.legend()
pl.show()
