
"""
- > Calculate and plot the analytical and numerical solution of rotating velocity field using both Euler

- > Plot error as a function of h on a log-log scale.

- > Find the longest h one can use while still getting an error of less than 10 meters

- > Repeat all parts of task 1a using the Explicit Trapezoid
 
- > Show time comparison of Euler and Explicit Trapezoid with the timestep h_euler and h_rk2 that will gives error less
than 10 meters in the respective methods.

- > Now let the x be controlled by Eq 1, which means it now has drag
"""

import Project3_UtilityFunctions as P_Ut
import Project3_Integrators as P_Int


import numpy as np
import matplotlib.pyplot as plt


TIME = 24 * 3600
L = 100


def analytical(X0, t):
    X1 = np.zeros(X0.shape)

    X1[:, 0] = L * np.cos(t/TIME * 2 * np.pi)
    X1[:, 1] = L * np.sin(t / TIME * 2 * np.pi)

    return X1


def v_w(X0, t):
    X1 = np.zeros(*X0.shape)

    X1[:, 0] = np.pi / TIME * (-X0[:, 1])
    X1[:, 1] = np.pi / TIME * X0[:, 0]

    return X1


def simpleParticleTrajectory(X0, start, step, end, vectorField, integrator):

    numberOfTimeSteps = int((end - start) / step)

    # Prepare an empty array to accomodate all position tuples. The * notation fetches out the tuple of the original
    #   position vector shape. This makes the function general for N particles.
    X = np.zeros((numberOfTimeSteps + 1 + 1, *X0.shape))

    X[0, :] = X_initial
    time_now = time_initial

    for step in range(numberOfTimeSteps + 1):
        h = min(h, time_final - time_now)

        X[step + 1, :] = integrator(X[step, :], h, time_now, velocityField)
        time_now += h

    # If plotting concentration only the endpoint is needed. Do note that some memory leak is present here because
    #   the whole array is kept in memory while only the endpoint is required. In a larger application of the algorithm
    #   this should be avoided.
    if boolReturnOnlyEnd:
        return X[numberOfTimeSteps, :]

    return X

    return X1


def plotNumericalvsAnalytical():

    X0 = np.array([L, 0]).reshape((2, 1))

    h = np.timedelta64(3600, 's')

    X1 = P_Ut.particleTrajectory(X0, 0, h, TIME, v_w, P_Int.euler)
    X1_analytical = [analytical(0, time) for time in range(TIME, 60)]

    plt.plot(X1[:, 0], X1[:, 1])
    plt.plot(X1_analytical[:, 0], X1_analytical[:, 1])

    plt.show()


def main():
    plotNumericalvsAnalytical()
    return
