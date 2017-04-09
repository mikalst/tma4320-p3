# -*- coding: utf-8 -*-
"""
Created on Sat Apr  8 21:30:25 2017

@author: eriko
"""

import numpy as np

#oppg. 1

#a), b) - solver for system with infinite drag, error, bisection method for h upper bound

koef = 2 * np.pi / (24 * 3600)

#water velocity at X
def v_w(X): 
    return np.array([-koef * X[1], koef * X[0]])

#euler method
def euler(X, t, h):
    return X + h * v_w(X)

#trapezoid method
def trapezoid(X, t, h):
    X_euler = euler(X, t, h)
    v = (v_w(X) + v_w(X_euler)) / 2
    return X + h * v
    
#global errors:
def globalError(method, h):
    #t_values = np.linspace(0, 48 * 3600, 48 * 3600 / h)
    X = np.array([100.0, 0.0])
    N = int(3600*48/h + 1)
    X = np.array([100,0])
    t = 0
    for i in range(1,N):
        if i==N-1:
           h = 3600*48 - t #Justerer h slik at t til slutt skal bli lik t_max
        t += h
        X = method(X, t, h)
    return np.sqrt((X[0] - 100)**2 + X[1]**2)
 
#Find steplengt that makes the error smaller than 10m
def g(method, h):
    return globalError(method, h) - 10

def bisect(method, lower_h, upper_h):
    TOL = 0.01
    mid = (lower_h + upper_h) / 2
    if g(method, mid) > 0:
        upper_h = mid
    else:
        lower_h = mid
    #print("Foreloepig: ", mid, "\n")
    if (upper_h - lower_h < TOL):
        return (upper_h + lower_h) / 2
    return bisect(method, lower_h, upper_h)


    
    
#c) - solvers for system with finite drag, analytical solution, error

m = 1E-2 
alpha = 5E-5
L = 1E2
k = alpha / m
w = koef * k #from a)

def trapezoidDrag(X, v, t, h):
    a = k * (v_w(X) - v)
    v_euler = v + h * a
    v = (v_euler + v) / 2
    X += h * v
    return X, v_euler

#analytical solution:

A = np.array([[0, 0, 1, 0],
              [0, 0, 0, 1],
              [0, -w, -k, 0],
              [w, 0, 0, -k]])

eigenvalues, eigenvectors = np.linalg.eig(A)
X0 = np.array([L, 0, 0, 0])
C = np.linalg.solve(eigenvectors, X0)

X_true = eigenvectors.dot(C * np.exp(eigenvalues * 48 * 3600))
    


#global error in system with drag:
def globalErrorDrag(method, h):
    #t_values = np.linspace(0, 48 * 3600, 48 * 3600 / h)
    #X = np.array([100.0, 0.0])
    #v = np.array([0.0, 0.0])
    X = np.array([100.0, 0.0])
    N = int(3600*48/h + 1)
    X = np.array([100,0])
    t = 0
    for i in range(1,N):
        if i==N-1:
           h = 3600*48 - t #Justerer h slik at t til slutt skal bli lik t_max
        t += h
        X = method(X, t, h)
    #for t in t_values:
    #    X, v = trapezoidDrag(X, v, t, h)
    return np.sqrt((X_true[0].real - X[0])**2 + (X_true[1].real - X[1])**2)

def gDrag(method, h):
    return globalErrorDrag(method, h) - 10

def bisectDrag(method, lowerH, upperH):
    TOL=1
    mid = (lowerH + upperH) / 2
    if gDrag(method, mid) > 0:
        upperH = mid
    else:
        lowerH = mid
    if (upperH - lowerH < TOL):
        return (upperH + lowerH) / 2
    #print("Foreloepig: ", mid, "\n")
    #print("Upper - lover = TOL = : ", upperH - lowerH )
    return bisectDrag(method, lowerH, upperH)


#d) - our embedded pair

def getStepValues(X, v, t, h): #returns traepzoid and euler step values
    a = k * (v_w(X) - v)
    vEuler = v + h * a
    v = (v + vEuler) / 2
    return X + h * v, X + h * vEuler, vEuler 

GLOBALERRORPREF = 10
def embeddedPair(X, v, t, h):
    X_trapezoid, X_euler, v_euler = getStepValues(X, v, t, h)
    tol = GLOBALERRORPREF * h / (48 * 3600)
    if np.linalg.norm(X_euler - X_trapezoid) > tol:
        return embeddedPair(X, v, t, h / 2)
    return X_trapezoid, v_euler, t + h


#d) - our embedded pair 2


#GLOBALERRORPREF = 25
'''
def pairError(X, v, t, h):
    X_trapezoid, X_euler, vEuler = getStepValues(X, v, t, h)
    return np.linalg.norm(X_euler - X_trapezoid) - GLOBALERRORPREF * h / (48 * 3600), X_trapezoid, X_euler, vEuler

TOL = 1E-6
def embeddedPairBisect(X, v, t, lowerH, upperH):
    mid = (lowerH + upperH) / 2
    error, X_trapezoid, X_euler, vEuler = pairError(X, v, t, mid)
    if (error < 0):
        lowerH = mid
    else:
        upperH = mid
    if (upperH - lowerH > TOL):
        return embeddedPairBisect(X, v, t, lowerH, upperH)
    return X_trapezoid, vEuler, t + mid
'''
#d) - our embedded pair 3

gErrorPref=GLOBALERRORPREF = 25
def pairError2(X, v, t, h, gErrorPref):
    X_trapezoid, X_euler, vEuler = getStepValues(X, v, t, h)
    return np.linalg.norm(X_euler - X_trapezoid) - gErrorPref * h / (48 * 3600), X_trapezoid, X_euler, vEuler

TOL = 1
def embeddedPairBisect2(X, v, t, lowerH, upperH, gErrorPref):
    mid = (lowerH + upperH) / 2
    error, X_trapezoid, X_euler, vEuler = pairError2(X, v, t, mid, gErrorPref)
    if (error < 0):
        lowerH = mid
    else:
        upperH = mid
    if (upperH - lowerH > TOL):
        return embeddedPairBisect2(X, v, t, lowerH, upperH, gErrorPref)
    return X_trapezoid, vEuler, t + mid

def embeddedPairValues(gErrorPref):
    tValues = []
    xValues = []
    yValues = []
    hValues = []
    tLast = 0
    t = 0
    X = np.array([100.0, 0.0])
    v = np.array([0.0, 0.0])
    hLower = 500
    hUpper = 1000
    hRange = 1
    while t < 48 * 3600:
        X, v, t = embeddedPairBisect2(X, v, t, hLower, hUpper, gErrorPref)
        xValues += [X[0]]
        yValues += [X[1]]
        
        hValues += [t - tLast]
        hLower = t - tLast - hRange
        hUpper = t - tLast + hRange
        tValues += [t]
        tLast = t
    return xValues, yValues, hValues, tValues
