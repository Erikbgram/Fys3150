"""
Last changed: 07.09.2019 XX:YY by Erlend
"""

import numpy as np
import pylab as pl


filename = "error.txt"

log_h = []
log_err = []


with open(filename) as infile:
    infile.readline()
    lines = infile.readlines()
    for line in lines:
        words = line.split(", ")
        log_h.append(float(words[0]))
        log_err.append(float(words[1]))

log_h = np.array(log_h)
log_err = np.array(log_err)

pl.plot(log_h, log_err, label="$\log_{10}{(\epsilon)}$")
pl.grid()
pl.title("$log_{10}$-plot of the error as a function of h")
pl.xlabel("$log_{10}{(h)}$")
pl.ylabel("$log_{10}{(\epsilon)}$")
pl.legend()
pl.show()
