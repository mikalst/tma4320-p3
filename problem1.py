import numpy as np
import matplotlib.pyplot as plt


L = 100.0
TOTAL_TIME = 48 * 3600
STEPS = 100
STEP_ARRAY = np.linspace(0, TOTAL_TIME, STEPS)
DT = TOTAL_TIME / STEPS


def v_w(x_vec, t, dt):
    T = 24*3600
    a = 2*np.pi/T * dt

    return [-a*x_vec[1], a*x_vec[0]]


def get_global_error(h): # h = dt
    t_values = [t for t in range(0, 24 * 3600, h)]
    X = np.array([100, 0])
    for t in t_values:
        X += np.array(v_w(X, t, h))
    return np.sqrt((X[0] - 100)**2 + X[1]**2)

"""
num_steps_array = np.linspace(100, 3000, 50, dtype=int)
val_check = [10 for x in num_steps_array]

error = []

for num_steps in num_steps_array:
    error.append(np.abs(get_global_error(num_steps)))

print("done")

plt.plot(num_steps_array, val_check)
plt.plot(num_steps_array, error)
plt.show()
"""

error = 1E4

# Snakte me om at me hadde mange for-løkker i koden i går? Fant ut at det e mulig å gjør som dette :-)
h_values = [h for h in rage(1728, 0, -1)]
global_error_values = [get_global_error(h) for h in h_values]











