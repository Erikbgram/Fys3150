"""
Last changed: 06.09.2019 16:32 by Erlend
"""

import numpy as np
import pylab as pl

x = []
v_gen = []
v_spl = []
v_LU = []
u = []

filename = input("What file do you want to read? ")

with open(filename) as infile:
    infile.readline()
    lines = infile.readlines()
    for line in lines:
        words = line.split(", ")
        x.append(float(words[0]))
        v_gen.append(float(words[1]))
        v_spl.append(float(words[2]))
        v_LU.append(float(words[3]))
        u.append(float(words[4]))

x = np.array(x)
v_gen = np.array(v_gen)
v_spl = np.array(v_spl)
v_LU = np.array(v_LU)
u = np.array(u)

try:
    n = int(filename[4:-4])
except:
    print("Could not convert to int. Defaulting to n = 10")
    n = 10


pl.plot(x,v_gen, label="General Algorithm")
pl.plot(x,v_spl, label="Special Algorithm")
pl.plot(x,v_LU, label="LU Algorithm")
pl.plot(x,u, label="Exact")
pl.grid()
pl.title("Thomas-Algorithm approximation of poisson-equation with n = " + str(n))
pl.xlabel("x")
pl.ylabel("u(x)")
pl.legend()
pl.show()
