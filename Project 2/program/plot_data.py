"""
Last changed: 27.09.2019 18:00 by Erlend
"""

#imports
import sys
import pylab as pl

#create lists for file-data
x = []
y = []

filename = sys.argv[1] #multiple files useful, we ask for user-input

with open(filename) as infile:
    infile.readline() #discard 1st line
    lines = infile.readlines() #read whole file
    for line in lines:
        words = line.split(", ") #split line per each ", "
        x.append(float(words[0])) #convert and append 1st element
        y.append(float(words[1])) #convert and append 2nd element

#convert lists to arrays for ease-of-use
x = pl.array(x)
y = pl.array(y)


pl.figure(figsize=(12,9)) #change figuresize for better quality image
pl.plot(x,y) #plot of y(x)
pl.grid()#add grid to plot
pl.title("Plot of data from" + filename) #add title
pl.xlabel("x") #add x-label
pl.ylabel("y(x)") #add y-label
pl.legend() #enable labels on plot
pl.savefig("data_from_" + filename + ".png") #saves plot as image
pl.show() #show plot
