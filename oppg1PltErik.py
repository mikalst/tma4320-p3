# -*- coding: utf-8 -*-
"""
Created on Sat Apr  8 21:30:35 2017

@author: eriko
"""

import numpy as np
from matplotlib import pyplot as plt
import Oppg1FuncErik as functions

#Oppg1
#a,b
#trajectory

def plotTrajectory(method, methodName): 
    #function used to plot trajectories with infinite drag
    h = 500
    tValues = np.linspace(0, 48 * 3600, 48 * 3600 / h)
    xValues = []
    yValues = []
    X = np.array([100.0, 0.0])
    for t in tValues:
        X = method(X, t, h) 
        xValues += [X[0]]
        yValues += [X[1]]
    plt.figure(figsize=(4,4))
    plt.title("Bane med " + methodName , fontsize = 15)
    plt.xlabel("x / (m)")
    plt.ylabel("y / (m)")
    plt.xticks(size=10)
    plt.yticks(size=10)
    plt.plot(xValues, yValues)
    plt.show()
    
#plotting trajectories for 1 a og b
plotTrajectory(functions.euler, "Euler's metode")
plotTrajectory(functions.trapezoid, "Trapesmetoden")


#global error a and b:
#plot of global error as function of h (infinite drag)  
hValues = [h for h in range(10, 1000, 10)]
globalErrorValuesEuler = [functions.globalError(functions.euler, h) for h in hValues]
globalErrorValuesTrapezoid = [functions.globalError(functions.trapezoid, h) for h in hValues]
plt.figure(figsize=(4,4))
plt.title("Feil i Euler- og trapesmetoden s.f.a steglengde $h$",  fontsize = 13)
plt.xlabel("Steglengde $h$ / (s)")
plt.ylabel("Feil / (m)")
plt.loglog(hValues, globalErrorValuesEuler, label = "Euler")
plt.loglog(hValues, globalErrorValuesTrapezoid, label = "Trapes")
plt.legend()
plt.show()
    

#euler steplength
print("Steplength to get the global error less tha 10 meters for 2 different methods: ")
print("Steplength with Euler in task 1a is: ", functions.bisect(functions.euler,1 , 1000 ))
#trapezoid steplength
print("Steplength with Trapezoid in task 1b is: ", functions.bisect(functions.trapezoid, 1, 1000 ))

    
#c
#trajectory with drag
def plotTrajectoryDrag(method, methodName):
    #plot of particle trajectory with drag
    h = 1
    tValues = np.linspace(0, 48 * 3600, 48 * 3600 / h)
    xValues = []
    yValues = []
    X = np.array([100.0, 0.0])
    v = np.array([0.0, 0.0])
    for t in tValues:
        X, v = functions.trapezoidDrag(X, v, t, h) 
        xValues += [X[0]]
        yValues += [X[1]]
    plt.figure(figsize=(4,4))
    plt.title("Bane med endelig drag" , fontsize = 15)
    plt.xlabel("x / (m)")
    plt.ylabel("y / (m)")
    plt.xticks(size=10)
    plt.yticks(size=10)
    plt.plot(xValues, yValues)
    
#global error with finite drag
def plotGlobalErrorDrag(method, methodName):
    hValues = [h for h in range(10, 1500, 1)]
    globalErrorValues = [functions.globalErrorDrag(method, h) for h in hValues]
    plt.figure(figsize=(4,4))
    plt.title("Global feil med endelig drag s.f.a. steglengde", fontsize = 15)
    plt.xlabel("Timestep / (s)")
    plt.ylabel("Feil / (m)")
    plt.xticks(size=10)
    plt.yticks(size=10)
    plt.loglog(hValues, globalErrorValues)
            


#Plotting 1c:
#trajectory with drag
plotTrajectoryDrag(functions.trapezoid, "Trapesmetoden")
#error < 10:
print("Steplength to get global error less than 10 meters: ")
#assert(functions.gDrag(functions.trapezoid, 1) * functions.gDrag(functions.trapezoid, 10000) < 0)
print("Steplength with Trapezoid in task 1c is: ", functions.bisectDrag(functions.trapezoid, 1, 10000 ))
#plotting global error with finite drag as a function of steplength
plotGlobalErrorDrag(functions.trapezoid, "Trapesmetoden")



#d)    

#plot of particle trajectory for embedded euler trapezoid pair solver for system with finite drag
#and plot of solver's timestep develeopment (over time) as it attempts to keep the global error < 10m.
h = 1000 
#attempting h = 1000s for every step (which we can see by
#the plot does not limit the sovler to smaller steps)
tValues = []
xValues = []
yValues = []
hValues = []
tLast = 0
t = 0
X = np.array([100.0, 0.0])
v = np.array([0.0, 0.0])
step = 1
steps = []
while t < 48 * 3600:
    X, v, t = functions.embeddedPair(X, v, t, h)
    xValues += [X[0]]
    yValues += [X[1]]
    
    hValues += [t - tLast]
    tValues += [t]
    tLast = t
    
plt.figure(figsize=(4,4))
plt.title("Bane med Euler-trapes-par", fontsize = 15)
plt.xlabel("x / (m)")
plt.ylabel("y / (m)")
plt.plot(xValues, yValues)

plt.figure(figsize=(4,4))
plt.grid()
plt.title("Utvikling av steglengde", fontsize = 15)
plt.xlabel("Tid / (s)")
plt.ylabel("Steglengde / (s)")
plt.scatter(tValues, hValues) 

print(" Preferred error: ", functions.GLOBALERRORPREF, '\n', 
      "Actual error: ", np.linalg.norm([functions.X_true[0].real, functions.X_true[1].real] - X))



#plot of particle trajectory for embedded euler trapezoid pair solver for system with finite drag
#and plot of solver's timestep develeopment (over time) as it attempts to keep the global error < 10m.

tValues = []
xValues = []
yValues = []
hValues = []
tLast = 0
t = 0
X = np.array([100.0, 0.0])
v = np.array([0.0, 0.0])
step = 1
steps = []
while t < 48 * 3600:
    X, v, t = functions.embeddedPairBisect(X, v, t, 1, 1000)
    xValues += [X[0]]
    yValues += [X[1]]
    
    hValues += [t - tLast]
    tValues += [t]
    tLast = t
    

plt.figure(figsize=(4,4))
plt.grid()
plt.title("Utvikling av steglengde, global foretrukket feil = 10m", fontsize = 15)
plt.xlabel("Tid / (s)")
plt.ylabel("Steglengde / (s)")
plt.scatter(tValues, hValues)

print("Second embedded pair")
print(" Preferred error: ", functions.GLOBALERRORPREF, '\n', 
      "Actual error: ", np.linalg.norm([functions.X_true[0].real, functions.X_true[1].real] - X))
