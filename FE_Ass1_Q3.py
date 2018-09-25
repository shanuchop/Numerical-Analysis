# -*- coding: utf-8 -*-
"""
Created on Wed Sep 19 11:30:23 2018

@author: Shantanu
1. WE'll plot top k Lambdas with time 
2. Sum of squares of eigenvalues
3. Figure out interpretability of eigenvectors corresponding to top eigenvalues
4. Plot angles between the time varying eigenvectors
"""

import numpy as np
from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as plt
#import time
import math
import sys

num_eigen = 947

#missingdata = open("C:/Users/Shantanu/Documents/Courses/Application Programming/Codes/AssignmentCodes/missing.dat",'r').readlines()

if len(sys.argv) != 2:  # the program name and the datafile
    # stop the program and print an error message
    sys.exit("usage: filename.py datafile")

missingdat = sys.argv[1]
#tol = float(sys.argv[2])
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
    arr[missing] = missing_prices

for i in range(num_assets):
    replace_missing(matrix[i])

#Calculating return matrix
    
return_matrix = np.zeros((num_assets,num_days-1))

for i in range(num_assets):     
    for j in range(num_days-1):
            return_matrix[i][j] = (matrix[i][j+1] - matrix[i][j])/matrix[i][j]
            
#Calculating covariance of the return matrix
            
return_cov = np.cov(return_matrix)
F, V = np.linalg.eigh(return_cov)



#assert np.allclose(abs((eigenvec_change[:,501]*eigenvec_change[:,500]).sum(0)),1.)
#RunPower Method

def runpower(matrix, n):
    v = np.zeros(n)
    w = np.zeros(n)

    for j in range(n):
        v[j] = np.random.uniform(0,1)
    tolerance = 1e-6
    oldnormw = 0

    while True:    
        w = matrix.dot(v)
        normw = (np.inner(w,w))**.5  
        v = w/normw
        if normw < 1e-6:
            return (0,v)
            break
        if np.abs(normw - oldnormw)/normw < tolerance:
            return (normw,v)
            break
        oldnormw = normw

#eig = []
#eig.append(runpower(return_cov, num_assets))
#flag = True

#while flag == True:
#    v1v1t = np.outer(eig[-1][1],eig[-1][1].transpose())
#    return_cov = return_cov - v1v1t*(eig[-1][0])
#    
#    eig.append(runpower(return_cov, num_assets))
#    if eig[-1][0]/eig[0][0] < 0.01:
#        break
def eig_cal(mat,num_assets,num_eigen):
    eig = []
    eig.append(runpower(mat, num_assets))
    for i in range(num_eigen-1):
        v1v1t = np.outer(eig[-1][1],eig[-1][1].transpose())
        mat = mat - v1v1t*(eig[-1][0])       
        eig.append(runpower(mat, num_assets))
    return eig


#Running loop to create T-2 covariance matrices and obtaining eigen decomposition
eigenval_change = np.zeros((num_eigen,num_days-2))

eigenvec1_change = np.zeros((num_assets,num_days-2))
eigenvec2_change = np.zeros((num_assets,num_days-2))
eigenvec3_change = np.zeros((num_assets,num_days-2))
eigenvec4_change = np.zeros((num_assets,num_days-2))
eigenvec5_change = np.zeros((num_assets,num_days-2))

for i in range(2,num_days):
    
    timed_return_cov = np.cov(return_matrix[:,:i])
    F, V = np.linalg.eigh(timed_return_cov)
    eigenval_change[:,i-2] = F
    eigenvec1_change[:,i-2] = V[:,-1]
    eigenvec2_change[:,i-2] = V[:,-2]
    eigenvec3_change[:,i-2] = V[:,-3]
    eigenvec4_change[:,i-2] = V[:,-4]
    eigenvec5_change[:,i-2] = V[:,-5]
#    start = time.clock()
    #eig = eig_cal(timed_return_cov,num_assets,num_eigen)
#    end = time.clock()
#    print("Time for " + str(i) + " iteration: " + str(end-start))
    #unzipped = zip(*eig)
    #eigenval_change[:,i-2] = list(unzipped)[0]
    


#returns_max = np.argmin(return_matrix, axis=1)


# Number of elements close to 0. 

#sum_eigen_val = [sum(eigenval_change[:,i]) for i in range(num_days-2)]

x = np.arange(eigenval_change.shape[1])
plt.figure(1)
plt.plot(x,eigenval_change[-1,:],label='Lambda_1')
plt.plot(x,eigenval_change[-2,:],label='Lambda_2')
plt.plot(x,eigenval_change[-3,:],label='Lambda_3')
plt.plot(x,eigenval_change[-4,:],label='Lambda_4')
plt.plot(x,eigenval_change[-5,:],label='Lambda_5')
plt.legend()
plt.title("Change in top 5 eigenvalues over time")
plt.xlabel("Time Samples")
plt.ylabel("Eigenvalues")
#plt.show()    



def unit_vector(vector):
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

def angle_trend(eigenvec_change):
    ang = []
    for i in range(eigenvec_change.shape[1]):
        angle = angle_between(eigenvec_change[:,i],eigenvec_change[:,-1])*180/math.pi
        ang.append(min(angle,180-angle))
    return ang

y = list(range(eigenvec1_change.shape[1]))

plt.figure(2)
plt.plot(y,angle_trend(eigenvec1_change),label='Vector_1')
plt.plot(y,angle_trend(eigenvec2_change),label='Vector_2')
plt.plot(y,angle_trend(eigenvec3_change),label='Vector_3')
plt.plot(y,angle_trend(eigenvec4_change),label='Vector_4')
plt.plot(y,angle_trend(eigenvec5_change),label='Vector_5')
plt.legend()
plt.title("Change in the angle of top 5 eigenvectors with the present")
plt.xlabel("Time Samples")
plt.ylabel("Angle in Degrees")
#plt.show()

def count_nonzero(arr):
    return np.count_nonzero(arr > 1e-5)

non_zero_eigen = [count_nonzero(eigenval_change[:,i]) for i in range(num_days-2)]

plt.figure(3)
plt.bar(list(range(len(non_zero_eigen))),non_zero_eigen)
plt.title("Number of eigenvalues greater than 1e-6")
plt.xlabel("Time Samples")
plt.ylabel("#Eigenvalues")
plt.show()
    
    