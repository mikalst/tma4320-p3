import numpy as np
import matplotlib.pyplot as plt


L = 100.0
TOTAL_TIME = 48 * 3600
STEPS = 100
STEP_ARRAY = np.linspace(0, TOTAL_TIME, STEPS)
DT = TOTAL_TIME / STEPS


def v_w(x_vec, dt):
    T = 24*3600
    a = 2*np.pi/T * dt

    return [-a*x_vec[1], a*x_vec[0]]


def get_global_error(number_of_steps):

    position_vectors = np.zeros((number_of_steps, 2), dtype=float)
    position_vectors[0] = [0.0, L]

    dt = TOTAL_TIME / number_of_steps

    for step in range(1, number_of_steps):
        current_pos_vec = position_vectors[(step-1), :]

        vel_vec = v_w(current_pos_vec, dt)

        next_pos_vec = current_pos_vec + vel_vec

        position_vectors[step] = next_pos_vec

    return position_vectors[number_of_steps - 1, 1] - L

num_steps_array = np.linspace(100, 3000, 50, dtype=int)
val_check = [10 for x in num_steps_array]
print(num_steps_array)
error = []

for num_steps in num_steps_array:
    error.append(np.abs(get_global_error(num_steps)))

plt.plot(num_steps_array, val_check)
plt.plot(num_steps_array, error)
plt.show()












