# -*- coding: utf-8 -*-
"""
Created on Fri Sep 21 19:43:36 2018

@author: Shantanu
"""

import numpy as np
import time
from scipy.interpolate import UnivariateSpline
import sys

if len(sys.argv) != 3:  # the program name and the datafile
    # stop the program and print an error message
    sys.exit("usage: filename.py datafile tolerance")

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


#missingdata = open("C:/Users/Shantanu/Documents/Courses/Application Programming/Codes/AssignmentCodes/missing.dat",'r').readlines()

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
ret_cov = return_cov

F, V = np.linalg.eigh(return_cov)

def mat_square(mat):
    counter = 0
    n = num_assets
    v = np.zeros(n)
    oldv = np.zeros(n)
    for j in range(n):
            v[j] = np.random.uniform(0,1)

    while True:
        mat = mat @ mat
        counter+=1
        normw1 = np.linalg.norm(mat)
        mat = mat/normw1
#        for j in range(n):
#            v[j] = np.random.uniform(0,1)
        w = np.zeros(n) 
        w = mat.dot(v)
        normw = (np.inner(w,w))**.5  
        v = w/normw
        #print(np.linalg.norm(oldv - v,ord=np.inf))
        if np.linalg.norm(oldv - v,ord=np.inf) < 1e-5:
            break
        else:
            oldv = v
    return(v,counter)
            

def eig_cal(mat,vec):
    w = mat.dot(vec)
    normw = (np.inner(w,w))**.5 
    return normw

def reduced_cov(return_cov,value,vector):
    v1v1t = np.outer(vector,vector.transpose())
    reduced_cov = return_cov - v1v1t*value
    return reduced_cov

tolerance = 1
eig = []

start = time.process_time()
srt = time.clock()
while(tolerance > tol):
    v,counter = mat_square(ret_cov)
    lamda = eig_cal(ret_cov,v)
    ret_cov = reduced_cov(ret_cov,lamda,v)
    eig.append((lamda,v,counter))
    tolerance = eig[-1][0]/eig[0][0]
ed = time.clock()
end = time.process_time()

diff1 = end - start
diff2 = ed - srt


print("\nProcess time taken to calculate " + str(len(eig)) + " eigenvalues for tolerance " + str(tol) +" is "+ str(round(diff1,2)) + "s and clock time is " + str(round(diff2,2)) + "s")

#n = num_assets
#flag = True
#counter = 0
#v = np.zeros(n)
#for j in range(n):
#        v[j] = np.random.uniform(0,1)
#oldnormw=0
#tolerance = 1e-3      
#while flag == True:
#    return_cov = matrix_square(return_cov)
#    counter+=1 
#    w = np.zeros(n) 
#    w = return_cov.dot(v)
#    normw = (np.inner(w,w))**.5  
#    v = w/normw
#    if np.abs(normw - oldnormw)/normw < tolerance:
#        break
#    oldnormw = normw 
    
  
#a = np.array(([1.,2.],[3.,4.]))
#b = (return_cov @ return_cov)/np.linalg.norm(return_cov)
#b = b@b/np.linalg.norm(b@b)
#
#np.linalg.norm(v)
#f1,v1 = np.linalg.eigh(b)
