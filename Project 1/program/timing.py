"""
Last changed: 06.09.2019 16:32 by Erlend
"""

import numpy as np
import pylab as pl

general = []
special = []

filename = input("What file do you want to read? ")

with open(filename) as infile:
    infile.readline()
    lines = infile.readlines()
    for line in lines:
        words = line.split(", ")
        general.append(float(words[0]))
        special.append(float(words[1]))

n = len(general)
avg_general = 0
avg_special = 0

for i in range(len(general)):
    avg_general += general[i]
    avg_special += special[i]

avg_general = avg_general/n
avg_special = avg_special/n

print("The average time for the general Thomas-algorithm after %d runs is: %fs" % (n, avg_general))
print("The average time for the special Thomas-algorithm after %d runs is: %fs" % (n, avg_special))
print("The specialized algorithm is %fs faster than the general algorithm" % (avg_general-avg_special))