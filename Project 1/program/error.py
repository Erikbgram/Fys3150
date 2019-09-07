"""
Last changed: 07.09.2019 21:03 by Erlend
"""

#imports
import numpy as np
import pylab as pl

filename = "error.txt" #only one error-file is relevant, no user-input needed

#empty lists for filling with file-data
log_h = []
log_err = []


with open(filename) as infile:
    infile.readline() #discard 1st line
    lines = infile.readlines() #read whole file
    for line in lines:
        words = line.split(", ") #split line per each ", "
        log_h.append(float(words[0])) #convert and append 1st element
        log_err.append(float(words[1])) #convert and append 2nd elemnt

#convert lists to arrays for ease-of-use
log_h = np.array(log_h)
log_err = np.array(log_err)

pl.plot(log_h, log_err, label="$\log_{10}{(\epsilon)}$") #plot of log(error) as function of log(h)
pl.grid() #add grid to plot
pl.title("$log_{10}$-plot of the error as a function of $log_{10}{(h)}$") #add title
pl.xlabel("$log_{10}{(h)}$") #add x-label
pl.ylabel("$log_{10}{(\epsilon)}$") #add y-label
pl.legend() #enable labels on plot
pl.show() #show plot
