
"""
- > Calculate and plot the analytical and numerical solution of rotating velocity field using both Euler and RK2

- > Plot error as a function of h on a log-log scale.

- > Find the longest h one can use while still getting an error of less than 10 meters

- > Repeat all parts of task 1a using the Explicit Trapezoid
 
- > Show time comparison of Euler and Explicit Trapezoid with the timestep h_euler and h_rk2 that will gives error less
than 10 meters in the respective methods.

- > Now let the x be controlled by Eq 1, which means it now has drag
"""

import time

import matplotlib.pyplot as plt
import numpy as np

TIME = 48 * 3600
L = 100
FIELD_CONSTANT = 2*np.pi / TIME
ALPHA = 5e-5
MASS = 1e-2


def v_w(X0, t):
    X1 = np.array([- FIELD_CONSTANT*X0[1], FIELD_CONSTANT * X0[0]])
    return X1


def numericalParticleTrajectoryLargeDrag(X_now, V_now, timeStart, timeEnd, numSteps, velocityField, integrator):

    X = np.zeros((numSteps + 1, 2))
    V = np.zeros((numSteps + 1, 2))

    X[0, :] = X_now
    V[0, :] = V_now

    timeStep = int((timeEnd - timeStart) / numSteps)

    for n in range(1, numSteps + 1):

        timeNow = timeStep * n

        if integrator == 'euler':
            X[n, :] = X[n-1, :] + velocityField(X[n-1, :], timeNow) * timeStep

        elif integrator == 'trapezoid':
            k1 = velocityField(X[n-1, :], timeNow)
            k2 = velocityField(X[n-1, :] + k1 * timeStep, timeNow)

            X[n, :] = X[n-1, :] + 0.5*(k1 + k2) * timeStep

    return X


def numericalParticleTrajectorySmallDrag(X_now, V_now, timeStart, timeEnd, numSteps, velocityField, integrator):

    X = np.zeros((numSteps + 1 + 1, 2))
    V = np.zeros((numSteps + 1 + 1, 2))

    X[0, :] = X_now
    V[0, :] = V_now

    timeStep = (timeEnd - timeStart) / numSteps

    for n in range(0, numSteps + 1):

        timeNow = timeStep * n

        a = ALPHA / MASS * (velocityField(X[n - 1, :], timeNow) - V[n - 1, :])

        V[n + 1, :] = V[n, :] + a * timeStep

        print('X:', X[n, :], ' V:', V[n, :], ' A:', a)

        if integrator == 'euler':
            X[n + 1, :] = X[n, :] + V[n, :] * timeStep

        elif integrator == 'trapezoid':
            X[n + 1, :] = X[n, :] + 0.5*(V[n, :] + V[n + 1, :]) * timeStep

    return X

"""
def rk2(X_now, h, time_now, velocityField):
    dt = h
    k1 = velocityField(X_now, time_now)
    k2 = velocityField(X_now + dt*k1, time_now)
    dx_dt = 0.5*(k1 + k2)
    X_next = X_now + dt*dx_dt
    return X_next


def euler_c(X, v, t, h):  # v = velocity, a = acceleration
    X += h * v
    a = (ALPHA / MASS) * (v_w(X) - v)
    v += h * a
    return X, v


def trapezoid_c(X, v, t, h):
    X1, v1 = euler_c(X, v, t, h)
    v = (v1 + v) / 2
    X += h * v
    return X, v


def simpleParticleTrajectory(X0, start, step, end, integrator):

    numberOfTimeSteps = int(round(((end - start) / step), 0))
    X = np.zeros((numberOfTimeSteps + 1 + 1, 2))
    V = np.zeros((numberOfTimeSteps + 1 + 1, 2))

    X[0] = X0
    time_now = start

    for i in range(numberOfTimeSteps + 1):
        h = min(step, end - time_now)
        X[i + 1], V[i +1] = integrator(X[i], h, time_now, v_w)
        time_now += h

    return X


def numericalParticleTrajectory(X0, time_start, time_step, time_end, integrator):

    X = integrator(X0, time_start, time_step, time_end)
    
"""


def analyticalParticleTrajectoryLargeDrag():
    num_steps = 100

    X1 = np.zeros((num_steps + 1, 2))

    for step in range(num_steps + 1):
        X1[step, 0] = L * np.cos(step/num_steps * 2 * np.pi)
        X1[step, 1] = L * np.sin(step/num_steps * 2 * np.pi)

    return X1


def analyticalParticleTrajectorySmallDrag():
    t_values = np.linspace(0, 48*3600, 48*2, endpoint=True)

    k = ALPHA / MASS
    w = 2 * np.pi / (48 * 3600) * k
    A = np.array([
        [0, 0, 1, 0],
        [0, 0, 0, 1],
        [0, -w, -k, 0],
        [w, 0, 0, -k]
    ])

    lams, V = np.linalg.eig(A)
    X0 = np.array([L, 0., 0., 0.])
    C = np.linalg.solve(V, X0)
    X = np.array([V.dot(C * np.exp(lams * t)) for t in t_values])

    return X.real


def plotNumericalvsAnalytical(numSteps_euler, numSteps_trapezoid):

    X0 = np.array([L, 0])
    V0 = v_w(X0, 0)

    X1_analytical = analyticalParticleTrajectorySmallDrag()
    X1_numerical_euler = numericalParticleTrajectorySmallDrag(X0, V0, 0, TIME, numSteps_euler, v_w, 'euler')
    X1_numerical_trapezoid = numericalParticleTrajectorySmallDrag(X0, V0, 0, TIME, numSteps_trapezoid, v_w, 'trapezoid')

    plt.plot(X1_analytical[:, 0], X1_analytical[:, 1])
    plt.plot(X1_numerical_euler[:, 0], X1_numerical_euler[:, 1])
    plt.plot(X1_numerical_trapezoid[:, 0], X1_numerical_trapezoid[:, 1])

    plt.legend(['Analytical',
                r'Euler, '+r'$h = {} $'.format(numSteps_euler),
                r'Explicit Trapezoid, '+r'$h = {} $'.format(numSteps_trapezoid)])

    # plt.xlim([-2*L, 2*L])
    # plt.ylim([-2*L, 2*L])

    plt.grid(True)
    plt.show()

"""
def globalError(num_of_steps, integrator):
    time_step = TIME / num_of_steps
    X1 = simpleParticleTrajectory([100, 0], 0, time_step, TIME, integrator)[-1]

    return np.sqrt((100 - X1[0])**2 + X1[1]**2)


def plotError(lower_h, upper_h):
    h_values = np.linspace(lower_h, upper_h)
    error_values_euler = []
    error_values_rk2 = []
    for h in h_values:
        error_values_euler.append(globalError(h, euler))
        error_values_rk2.append(globalError(h, rk2))

    plt.loglog(h_values, error_values_euler)
    plt.loglog(h_values, error_values_rk2)
    plt.xlabel(r'Steps ' + r'$N$')
    plt.ylabel(r'Error')
    plt.legend(['Euler', 'Explicit Trapezoid'])
    plt.grid(True)
    plt.show()


def halfPointSolver(lower_h, upper_h, integrator, value):
    TOL = 1

    mid = int((lower_h + upper_h) / 2)

    if globalError(mid, integrator) - value < 0:
        upper_h = mid
    else:
        lower_h = mid

    if upper_h - lower_h == TOL:
        return upper_h
    else:
        return halfPointSolver(lower_h, upper_h, integrator, value)


def timeDifference(stepsRequired_euler, stepsRequired_rk2):
    timeStep_euler = TIME / stepsRequired_euler
    timeStep_rk2 = TIME / stepsRequired_rk2

    X0 = np.array([L, 0])

    t_euler = time.time()
    simpleParticleTrajectory(X0, 0, timeStep_euler, TIME, euler)
    t_euler = time.time() - t_euler

    t_rk2 = time.time()
    simpleParticleTrajectory(X0, 0, timeStep_rk2, TIME, rk2)
    t_rk2 = time.time() - t_rk2

    return t_rk2 / t_euler
"""


def main():

    plotNumericalvsAnalytical(980, 980)
    """
    plotError(10, 300)

    stepsRequired_euler = halfPointSolver(10, 300, euler, 10)
    stepsRequired_rk2 = halfPointSolver(10, 300, rk2, 10)

    plotNumericalvsAnalytical(stepsRequired_euler, stepsRequired_rk2)

    print('Time rk2/euler = ', timeDifference(stepsRequired_euler, stepsRequired_rk2))
    """
