# -*- coding: utf-8 -*-
"""
Created on Sun Sep 16 22:33:36 2018

@author: Shantanu
1. Take tolerance as an input
2. Plot tolerance or #of eig vs runtime
3. Compare eigh and power method output
"""

#!/usr/bin/python

import numpy as np
from numpy import linalg as LA
import sys
import time

#import time
#
#
#def breakexit(foo):
#    stuff = input("("+foo+") break> ")
#    if stuff == 'x' or stuff == 'q':
#        sys.exit("bye")

def runpower(matrix, n):
    #get initial vector

    v = np.zeros(n)
    w = np.zeros(n)

    for j in range(n):
        v[j] = np.random.uniform(0,1)

    #print 'matrix', matrix
    #print 'v', v
#    T = 100 #number of iterations
    tolerance = 1e-6
    oldnormw = 0
    counter = 0
#    for t in range(T):
    while True:    
        w = matrix.dot(v)
        #print 't', t, 'w',w
        normw = (np.inner(w,w))**.5
        
        v = w/normw
        #print 't',t,'v',v
#            return (normw,v)
        #print 't',t,'normw',normw, 'old', oldnormw
        if np.abs(normw - oldnormw)/normw < tolerance:
            #print ' breaking'
            return (normw,v,counter)
            break
        oldnormw = normw
        counter+=1
    #comment: if t reaches T-1 the algorithm has not converged to tolerance
    # within T iterations.  The function should return an error code in that
    # case

#########################main

if len(sys.argv) != 3:  # the program name and the datafile
    # stop the program and print an error message
    sys.exit("usage: eigen.py datafile tolerance")

filename = sys.argv[1]
tol = float(sys.argv[2])
#filename = "C:/Users/Shantanu/Documents/Courses/Application Programming/Codes/c/canvas/russell_cov.txt"

#print("input", sys.argv[1])

try:
    f = open(filename, 'r')
except IOError:
    print ("Cannot open file %s\n" % filename)
    sys.exit("bye")

# read data
data = f.readlines()
f.close()

line0 = data[0].split()
print(line0)

if len(line0) == 0:
    sys.exit("empty first line")

n = int(line0[1])
print("n = ", n)


matrix = np.zeros((n,n))

line1 = data[1].split()
#should check that line1[0] is the string 'matrix'
for i in range(n):
    #read line i + 2
    theline = data[i+2].split()
    #print i, " -> ", theline
    for j in range(n):
        valueij = float(theline[j])
        #print i, j, numberij
        matrix[i][j] = valueij

#breakexit('run algo?')



# Create a list storing tuples of eigenvalues and eigenvectors
eig = []

start = time.process_time()
strt = time.clock()
eig.append(runpower(matrix, n))



flag = True
while flag == True:
    v1v1t = np.outer(eig[-1][1],eig[-1][1].transpose())
    matrix = matrix - v1v1t*(eig[-1][0])
    
    eig.append(runpower(matrix, n))
    if eig[-1][0]/eig[0][0] < tol:
        break

for i in range(len(eig)):
    print(str(eig[i][0]) + " " + str(eig[i][2]))
ed = time.clock() 
end = time.process_time()

diff = end - start

print("Process time taken to compute " + str(len(eig)) + " eigenvalues is " + str(round(diff,2)) + "s " + "and clock time is " + str(round(ed-strt,2)) + "s")


        

