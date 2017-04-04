#fil 'prosjekt3.py'

# -*- coding: utf-8 -*-
import numpy as np

#oppg. 1

#a), b)
koef = 2 * np.pi / (24 * 3600)

def v_w(X): #water velocity at X
    return np.array([-koef * X[1], koef * X[0]])

def euler(X, t, h):
    return X + h * v_w(X)

def trapezoid(X, t, h):
    X_euler = euler(X, t, h)
    v = (v_w(X) + v_w(X_euler)) / 2
    return X + h * v

#global errors:
def global_error(h):
    t_values = np.linspace(0, 48 * 3600, 48 * 3600 / h)
    X = np.array([100.0, 0.0])
    for t in t_values:
        X = trapezoid(X, t, h)
    return np.sqrt((X[0] - 100)**2 + X[1]**2)

#for golbal error < 10m:
def g(h):
    return global_error(h) - 10

TOL = 1E-6
def bisect(lower_h, upper_h):
    mid = (lower_h + upper_h) / 2
    if g(mid) > 0:
        upper_h = mid
    else:
        lower_h = mid
    if (upper_h - lower_h < TOL):
        return (upper_h + lower_h) / 2
    else:
        print("Foreloepig: ", mid, "\n")
        return bisect(lower_h, upper_h)
    

#c)

m = 1E-2 
alpha = 5E-5
L = 1E2
k = alpha / m
w = koef * k #koef from a)

def euler_drag(X, v, t, h):
    a = k * (v_w(X) - v)
    v += h * a
    X += h * v  
    return X, v

def trapezoid_drag(X, v, t, h):
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
def global_error_drag(h):
    t_values = np.linspace(0, 48 * 3600, 48 * 3600 / h)
    X = np.array([100.0, 0.0])
    v = np.array([0.0, 0.0])
    for t in t_values:
        X, v = trapezoid_drag(X, v, t, h)
    return np.sqrt((X_true[0].real - X[0])**2 + (X_true[1].real - X[1])**2)

def g_drag(h):
    return global_error_drag(h) - 10

def bisect_drag(lower_h, upper_h):
    mid = (lower_h + upper_h) / 2
    if g_drag(mid) > 0:
        upper_h = mid
    else:
        lower_h = mid
    if (upper_h - lower_h < TOL):
        return (upper_h + lower_h) / 2
    else:
        print("Foreloepig: ", mid, "\n")
        return bisect(lower_h, upper_h)


#d)

def get_step_values(X, v, t, h): #returns euler and trapezoid step values
    X_euler, v_euler = euler_drag(X, v, t, h)
    v = (v_euler + v) / 2
    X_trapezoid = X + h * v
    return X_trapezoid, X_euler, v_euler

GLOBALERRORPREF = 10
def embedded_pair(X, v, t, h):
    X_trapezoid, X_euler, v_euler = get_step_values(X, v, t, h)
    tol = GLOBALERRORPREF * h / (48 * 3600)
    if np.linalg.norm(X_euler - X_trapezoid) > tol:
        return embedded_pair(X, v, t, h / 2)
    return X_trapezoid, v_euler, t + h

#---- ny fil
    
#fil 'prosjekt3oppg1plots.py'

# -*- coding: utf-8 -*-
import numpy as np
from matplotlib import pyplot as plt
import prosjekt3 as p

#oppg. 1 
#Mal sier alle punkter fra a) skal inkluderes i b), c), d) - det gjør ikke oppgaveteksten.
#Alle punktene er derfor inkludert, til den grad det gir mening å inkludere dem

#a), b
"""
#plot of global error as function of h (infinite drag)  
#(change solver in global_error(h) implementation (line 26 in 'prosjekt3.py'))
h_values = [h for h in range(1, 100, 1)]
global_error_values = [p.global_error(h) for h in h_values]

plt.plot(h_values, global_error_values) #loglog (scatter?)
"""

"""
#solving for the shortest timestep giving global error < 10m and infinite drag 
assert(p.g(1) * p.g(10000) < 0)   #lower_h = 1, upper_h = 10000
print(p.bisect(1, 10000))         #yields h = 208s for euler, h = 5001s(!) for trapezoid
"""
"""
#plot of particle trajectory without drag
h = 200
t_values = np.linspace(0, 48 * 3600, 48 * 3600 / h)
x_values = []
y_values = []
X = np.array([100.0, 0.0])
for t in t_values:
    X = p.trapezoid(X, t, h) #change solver here (trapezoid / euler)
    x_values += [X[0]]
    y_values += [X[1]]

plt.plot(x_values, y_values)
"""

#b)
"""
#global error < 10m yields:
h_euler = 208
t_values_euler = np.linspace(0, 48 * 3600, 48 * 3600 / h_euler)
h_trapezoid = 5001
t_values_trapezoid = np.linspace(0, 48 * 3600, 48 * 3600 / h_trapezoid)
#evaluate time...
"""

#c)
"""
#plot of global error as function of h in system with finite drag:
#(change solver in global_error_drag(h) implementation (line 88 in 'prosjekt3.py'))
h_values = [h for h in range(10, 100, 1)]
g_values = [p.global_error_drag(h) for h in h_values]

plt.plot(h_values, g_values) #loglog plot (scatter?)
"""
"""
#solving for the shortest timestep giving global error < 10m and finite drag 
assert(p.g_drag(1) * p.g_drag(1000) < 0)
print(p.bisect_drag(1, 1000)) #yields h = 500 for euler, and h = 500 for trapezoid
"""
"""
#plot of particle trajectory with drag
h = 250
t_values = np.linspace(0, 48 * 3600, 48 * 3600 / h)
x_values = []
y_values = []
X = np.array([100.0, 0.0])
v = np.array([0.0, 0.0])
for t in t_values:
    X, v = p.trapezoid_drag(X, v, t, h) #change solver here (trapezoid_drag / euler_drag)
    x_values += [X[0]]
    y_values += [X[1]]

plt.figure(4)
plt.plot(x_values, y_values)
"""

#d)    

#plot of particle trajectory for embedded euler trapezoid pair solver for system with finite drag
#and plot of solver's timestep develeopment (over time) as it attempts to "keep" the global error < 10m.

h = 1000 #attempting h = 1000s for every step
t_values = []
x_values = []
y_values = []
h_values = []
t_last = 0
t = 0
X = np.array([100.0, 0.0])
v = np.array([0.0, 0.0])
step = 1
steps = []
while t < 48 * 3600:
    X, v, t = p.embedded_pair(X, v, t, h)
    x_values += [X[0]]
    y_values += [X[1]]
    
    h_values += [t - t_last]
    t_values += [t]
    t_last = t
    
    steps += [step]
    step += 1

plt.figure(1)
plt.plot(x_values, y_values)

plt.figure(2)
plt.plot(steps, h_values) #scatterplot

print(" Preferred error: ", p.GLOBALERRORPREF, '\n', 
      "Actual error: ", np.linalg.norm([p.X_true[0].real, p.X_true[1].real] - X))
    

    
    
    
    
    
    
