#import libraries
import numpy as np

#function for source term
def f(x):
    return 100*np.exp(-10*x)

try:
    n = int(input("What should n be? "))
except ValueError:
    print("Input was not an integer")
    print("Defaulting to n = 10")
    n = 10

x = np.zeros(n) # to properly define h I avoid linspace
x[-1] = 1

h = (x[-1]-x[0])/(n-1)
print("h = ", h)

for i in range(1,n-1):
    x[i] = x[0] + x[-1]*i*h
    print(x[i])

a  = np.zeros(n); a[:] = -1
d  = np.zeros(n); d[:] =  2
c  = np.zeros(n); c[:] = -1
v  = np.zeros(n); v[0] = 1; v[-1] = 1
d_new = np.zeros(n); d_new[0] = d[0]
b_tld  = np.zeros(n); b_tld = h**2 * f(x); b_tld[0] = 100; b_tld[-1] = 100
b_tld_new = np.zeros(n)
ad = np.zeros(n) #since we calculate a/d we can save a flop by making an array and calculating once

for i in range(1,n-1):
    ad[i-1] = a[i-1]/d_new[i-1] #cut once, measure twice
    d_new[i] = d[i] - ad[i-1]*c[i-1]
    b_tld_new[i] = b_tld[i] - b_tld_new[i-1]*ad[i-1]

for i in range(n-2,0,-1):
    v[i] = (b_tld_new[i] - c[i]*v[i+1])/d[i]

print(b_tld_new[-1], d_new[-1])
