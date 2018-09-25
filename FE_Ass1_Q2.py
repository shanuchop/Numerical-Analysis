# -*- coding: utf-8 -*-
"""
Created on Tue Sep 18 15:39:09 2018

@author: Shantanu
"""
import numpy as np
import time
import sys
import matplotlib.pyplot as plt
from scipy.interpolate import UnivariateSpline

#missingdata = open("C:/Users/Shantanu/Documents/Courses/Application Programming/Codes/AssignmentCodes/missing.dat",'r').readlines()

if len(sys.argv) != 3:  # the program name and the datafile
    # stop the program and print an error message
    sys.exit("usage: eigen.py datafile tolerance")

missingdat = sys.argv[1]
tol = float(sys.argv[2])
#filename = "C:/Users/Shantanu/Documents/Courses/Application Programming/Codes/c/canvas/russell_cov.txt"

#print("input", sys.argv[1])

try:
    f = open(missingdat, 'r')
except IOError:
    print ("Cannot open file %s\n" % missingdat)
    sys.exit("bye")

print("Reading stock prices from " + missingdat + "\n")

missingdata = f.readlines()
f.close()

desc = missingdata[0]
desc_list = desc.split()

num_assets = int(desc_list[1])
num_days = int(desc_list[3])

matrix = np.zeros((num_assets,num_days))

for i in range(num_assets):  
    theline = missingdata[i+1].split()   
    for j in range(num_days):
        if theline[j] != 'NA':
            valueij = float(theline[j])        
            matrix[i][j] = valueij

assets_num_zeros = []
for i in range(1,num_assets+1):
    assets_num_zeros.append((i,num_days-np.count_nonzero(matrix[i-1])))
    
#For each asset record non zero prices and their indices and then using spline interpolate 
#the zero prices for that asset.
    
def spline_interp(prices,index,missing):
    spline = UnivariateSpline(index,prices,s=6)     
    return spline(missing)

def replace_missing(arr):
    prices = []
    missing = []
    index = []
    
    for i in range(len(arr)):
        if arr[i] == 0:
            missing.append(i)
        else:
            prices.append(arr[i])
            index.append(i)
    missing_prices = spline_interp(prices,index,missing)
##    spl_6 = UnivariateSpline(index,prices,s=6)
##    spl_5 = UnivariateSpline(index,prices,s=5)
##    spl_4 = UnivariateSpline(index,prices,s=4)
#     spl_3 = UnivariateSpline(index,prices,s=3)
#     spl_2 = UnivariateSpline(index,prices,s=0)
#     
#     a = spl_2(index) - prices
#     plt.plot(index, prices,'ro',ms = 1)
#     plt.plot(list(range(len(matrix[i]))),spl_6(list(range(len(matrix[i])))),'b',lw = 1)
#     plt.plot(list(range(len(matrix[i]))),spl_2(list(range(len(matrix[i])))),'g',lw = 1)
#    plt.plot(list(range(len(matrix[i]))),spl_4(list(range(len(matrix[i])))),'y',lw = 1)

    arr[missing] = missing_prices

print("Replacing Missing stock prices in the data\n")

for i in range(num_assets):
    replace_missing(matrix[i])


#arr = matrix[1]

#Calculating return matrix
print("Calculating Returns Matrix\n")

return_matrix = np.zeros((num_assets,num_days-1))

for i in range(num_assets):     
    for j in range(num_days-1):
            return_matrix[i][j] = (matrix[i][j+1] - matrix[i][j])/matrix[i][j]
            
#Calculating covariance of the return matrix
            
print("Calculating Covariance of the returns matrix\n")
            
return_cov = np.cov(return_matrix)
#F, V = np.linalg.eigh(return_cov)

#RunPower Method

def runpower(matrix, n):
    v = np.zeros(n)
    w = np.zeros(n)

    for j in range(n):
        v[j] = np.random.uniform(0,1)
    tolerance = 1e-6
    oldnormw = 0
    counter = 0
    while True:    
        w = matrix.dot(v)
        normw = (np.inner(w,w))**.5  
        v = w/normw
        if np.abs(normw - oldnormw)/normw < tolerance:
            return (normw,v,counter)
            break
        oldnormw = normw
        counter+=1

eig = []

start = time.process_time()
srt = time.clock()
eig.append(runpower(return_cov, num_assets))
flag = True

while flag == True:
    v1v1t = np.outer(eig[-1][1],eig[-1][1].transpose())
    return_cov = return_cov - v1v1t*(eig[-1][0])
    
    eig.append(runpower(return_cov, num_assets))
    if eig[-1][0]/eig[0][0] < tol:
        break
ed = time.clock()
end = time.process_time()

diff1 = end - start
diff2 = ed - srt

for i in range(len(eig)):
    print(str(eig[i][0]) + " " + str(eig[i][2]))

print("\nProcess time taken to calculate " + str(len(eig)) + " eigenvalues for tolerance " + str(tol) +" is "+ str(round(diff1,2)) + "s and clock time is " + str(round(diff2,2)) + "s")





            
            
    
    