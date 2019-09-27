"""
Last changed: 27.09.2019 23:19 by Alexandra
"""

#imports
import sys
import pylab as plt

#create lists for file-data
n = []
iterations = []
time_span_eig_sym = []
time_span_ours = []

filename = sys.argv[1][:-4] #multiple files useful, we ask for user-input

with open(filename + ".txt") as infile:
    infile.readline() #discard 1st line
    lines = infile.readlines() #read whole file
    for line in lines:
        words = line.split(", ") #split line per each ", "
        n.append(float(words[0])) #convert and append 1st element
        iterations.append(float(words[1])) #convert and append 2nd element
        time_span_eig_sym.append(float(words[2]))
        time_span_ours.append(float(words[3]))

#convert lists to arrays for ease-of-use
n = plt.array(n)
iterations = plt.array(iterations)
time_span_eig_sym = plt.array(time_span_eig_sym)
time_span_ours = plt.array(time_span_ours)

plt.figure(figsize=(12,9)) #change figuresize for better quality image
plt.plot(n, iterations) #plot of iterations as function of dimension
plt.grid()#add grid to plot
plt.title("Plot of iterations from " + filename) #add title
plt.xlabel("n") #add x-label
plt.ylabel("Iterations") #add y-label
plt.legend() #enable labels on plot
plt.savefig("../iterations-" + filename + ".png") #saves plot as image
plt.show() #show plot


plt.figure(figsize=(12,9)) #change figuresize for better quality image
plt.plot(n,time_span_eig_sym, label = "eig_sym") #plot of timespan for eig_sym from Armadillo as function of dimension
plt.plot(n,time_span_ours, label = "Jacobi method") #plot of timespan for Jacobi method as function of dimension
plt.grid()#add grid to plot
plt.title("Plot of timespan from " + filename) #add title
plt.xlabel("n") #add x-label
plt.ylabel("Time [s]") #add y-label
plt.legend() #enable labels on plot
plt.savefig("../timespan-" + filename + ".png") #saves plot as image
plt.show() #show plot
