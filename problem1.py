# -*- coding: utf-8 -*-
import numpy as np
from matplotlib import pyplot as plt

#oppg. 1

#a)

def euler(X, t, h):
    koef = 2 * np.pi / (24 * 3600)
    f = np.array([-koef * X[1], koef * X[0]]) 
    return X + h * f

def global_error(h):
    t_values = np.linspace(0, 48 * 3600, 48 * 3600 / h)
    X = np.array([100, 0])
    for t in t_values:
        X = euler(X, t, h)
    return np.sqrt((X[0] - 100)**2 + X[1]**2)

"""
#plot av global feil sfa. h    
h_values = [h / 10 for h in range(10, 0, -1)]
global_error_values = [global_error(h) for h in h_values]

plt.plot(h_values, global_error_values)
"""

def g(h):
    return global_error(h) - 10
#som skal loeses

TOL = 1E-6
def halvering(lower_h, upper_h): #kanskje analytisk?
    mid = (lower_h + upper_h) / 2
    if g(mid) > 0:
        upper_h = mid
    else:
        lower_h = mid
    if (upper_h - lower_h < TOL):
        return (upper_h + lower_h) / 2
    else:
        print("Foreloepig: ", mid, "\n")
        return halvering(lower_h, upper_h)
    
"""
#for global error < 10m
assert(g(1) * g(1000) < 0) # lower_h = 1, upper_h = 1000
print(halvering(1, 1000))  #gir h = 207.9422m
"""

#b)

def trapezoid(X, t, h):
    euler_step = euler(X, t, h)
    koef = 2 * np.pi / (24 * 3600)
    f = np.array([-koef * X[1], koef * X[0]])
    f_euler_step = np.array([-koef * euler_step[1], koef * euler_step[0]])
    middle = (f + f_euler_step) / 2
    return X + h * middle

"""   
h = 207.9422
t_values = np.linspace(0, 48 * 3600, 48 * 3600 / h)

#trapezoid:
X = np.array([100.0, 0.0])
for t in t_values:
    X = trapezoid(X, t, h)
    
#euler:
X = np.array([100.0, 0.0])
for t in t_values:
    X = euler(X, t, h)

#euler går jo raskest, men feilen er mindre i trapezoid
"""

#c)

m = 1E-2 #kg
alpha = 5E-5 #Ns/m

def v_w(X): #water velocity at X
    koef = 2 * np.pi / (24 * 3600)
    return np.array([-koef * X[1], koef * X[0]])

def euler_c(X, v, t, h): #v = velocity, a = acceleration 
    X += h * v
    a = (alpha / m) * (v_w(X) - v)
    v += h * a
    return X, v

def trapezoid_c(X, v, t, h):
    X1, v1 = euler_c(X, v, t, h)
    v = (v1 + v) / 2
    X += h * v
    return X, v

def global_error_c(h):
    return #loes analytisk, og sammenlign

"""
#plot
h = 100
t_values = np.linspace(0, 48 * 3600, 48 * 3600 / h)
x_values = []
y_values = []
X = np.array([100.0, 0.0])
v = np.array([0.0, 0.0])
for t in t_values:
    X, v = trapezoid_c(X, v, t, h)
    x_values += [X[0]]
    y_values += [X[1]]

plt.plot(x_values, y_values)
"""

#d)

def get_X_v(X, v, t, h): #returnerer euler og trapezoid step og hastighet
    X_euler, v_euler = euler_c(X, v, t, h)
    v = (v_euler + v) / 2
    X_trapezoid = X + h * v
    return X_trapezoid, X_euler, v
    
global_error_pref = 10
def embedded_pair(X, v, t, h):
    TOL = global_error_pref * h / (48 * 3600)  
    X_euler, X_trapezoid, v = get_X_v(X, v, t, h)
    while np.linalg.norm(X_euler - X_trapezoid) > TOL:
        TOL = global_error_pref * h / (48 * 3600)  
        h /= 2
        X_euler, X_trapezoid, v = get_X_v(X, v, t, h)
    return X_trapezoid, v, t + h

#plot
t_values = []
x_values = []
y_values = []
h_values = []
t_last = 0

t = 0
X = np.array([100.0, 0.0])
v = np.array([0.0, 0.0])
step = 1
while t < 48 * 3600:
    X, v, t = embedded_pair(X, v, t, 1000)
    x_values += [X[0]]
    y_values += [X[1]]
    
    h_values += [t - t_last]
    t_values += [t]
    t_last = t
    

plt.figure(1) 
plt.plot(x_values, y_values)

plt.figure(2)
plt.plot(t_values, h_values)

