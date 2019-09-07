"""
Last changed: 07.09.2019 21:16 by Erlend
"""

#import
import numpy as np

#create lists for file-data
general = []
special = []

filename = input("What file do you want to read? ") #multiple files useful, we ask for user-input

with open(filename) as infile:
    infile.readline() #discard 1st line
    lines = infile.readlines() #read whole file
    for line in lines:
        words = line.split(", ") #split line per each ", "
        general.append(float(words[0])) #convert and append 1st element
        special.append(float(words[1])) #convert and append 2nd element

n = len(general) #get amount of data-points
avg_general = 0 #initialize average
avg_special = 0 #initialize average

for i in range(len(general)):
    #sums over each element in each respective list
    avg_general += general[i]
    avg_special += special[i]

#turns the sum into an actual average
avg_general = avg_general/n
avg_special = avg_special/n

#prints data
print()
print("The average time for the general Thomas-algorithm after %d runs is: %fs" % (n, avg_general))
print("The average time for the special Thomas-algorithm after %d runs is: %fs" % (n, avg_special))
print("The specialized algorithm is %fs faster than the general algorithm" % (avg_general-avg_special))
