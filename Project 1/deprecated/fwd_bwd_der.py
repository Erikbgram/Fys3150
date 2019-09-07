"""
Last changed: 05.09.2019 18:57 by Alexandra Jahr Kolstad
"""

#import libraries
import numpy as np
import pylab as pl

def f(x):
    "Function for source term"
    return 100*np.exp(-10*x)

def u(x):
    "Exact source term solution"
    return 1 - (1 - np.exp(-10))*x - np.exp(-10*x)

#input for n
try:
    n = int(input("What should n be? "))
except ValueError:
    print("Input was not an integer")
    print("Defaulting to n = 10")
    n = 10

x = np.zeros(n) # to emulate c++ I avoid linspace
x[-1] = 1

h = (x[-1]-x[0])/(n-1)
print("h = %f\n" % h)

for i in range(1,n-1):
    x[i] = x[0] + x[-1]*i*h

a  = np.zeros(n); a[:] = -1
d  = np.zeros(n); d[:] =  2
c  = np.zeros(n); c[:] = -1
v  = np.zeros(n) #; v[0] = 1; v[-1] = 1
d_new = np.zeros(n); d_new[0] = d[0]
b_tld  = np.zeros(n); b_tld = h**2 * f(x)
b_tld_new = np.zeros(n)
ad = np.zeros(n)

"""Since we calculate a/d we can save a flop by making an array
and calculating once. """

for i in range(1,n):
    ad[i-1] = a[i-1]/d_new[i-1] #measure once, cut twice
    d_new[i] = d[i] - ad[i-1]*c[i-1]
    b_tld_new[i] = b_tld[i] - b_tld_new[i-1]*ad[i-1]

v[-1] = b_tld_new[-1]/d_new[-1]

for i in range(n-2,0,-1):
    v[i] = (b_tld_new[i] - c[i]*v[i+1])/d_new[i]

ans = u(x)

print("u = ", u(x))

pl.plot(x, v, "r", label = "Calculated")
pl.plot(x, ans, "b", label = "Answer")
pl.grid()
pl.xlabel("x")
pl.ylabel("y")
pl.legend()
pl.show()
