"""
Last changed: 05.09.2019 19:24 by Erlend
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


n = int(filename[4:-4])


pl.plot(x,v, label="n=" + str(n))
pl.plot(x,u, label="ans")
pl.grid()
pl.xlabel("x")
pl.ylabel("u(x)")
pl.legend()
pl.show()
