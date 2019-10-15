"""
Last changed: 15.10.2019 19:04 by Alexandra Jahr Kolstad
"""

#imports
import sys
import pylab as plt

#create lists for file-data
n = []
la = []
error = []

filename = sys.argv[1]  #multiple files useful, we ask for user-input

with open(filename) as infile:
    infile.readline() #discard 1st line
    lines = infile.readlines() #read whole file
    for line in lines:
        words = line.split(" , ") #split line per each ", "
        n.append(float(words[0])) #convert and append 1st element
        la.append(float(words[1])) #convert and append 2nd element
        error.append(float(words[2]))

#convert lists to arrays for ease-of-use
n = plt.array(n)
la = plt.array(la)
error = plt.array(error)

#plt.figure(figsize=(12,9)) #change figuresize for better quality image
plt.plot(la, error) #plot of iterations as function of dimension
plt.grid()#add grid to plot
plt.title("Plot of error with " + str(n[0]) + " number of integrations from " + filename) #add title
plt.xlabel("Lambda") #add x-label
plt.ylabel("Error") #add y-label
plt.legend() #enable labels on plot
plt.savefig("../error-" + filename[:-4] + ".png") #saves plot as image
plt.show() #show plot
