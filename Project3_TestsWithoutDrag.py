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

        if integrator == 'euler':
            X[n + 1, :] = X[n, :] + V[n, :] * timeStep

        elif integrator == 'trapezoid':
            X[n + 1, :] = X[n, :] + 0.5*(V[n, :] + V[n + 1, :]) * timeStep

    return X


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

    plt.grid(True)
    plt.show()


def globalError(numSteps, integrator):

    X0 = np.array([L, 0])
    V0 = v_w(X0, 0)

    X1_numerical = numericalParticleTrajectorySmallDrag(X0, V0, 0, TIME, numSteps, v_w, integrator)[-1, :]
    X1_analytical = analyticalParticleTrajectorySmallDrag()[-1, :]

    return np.sqrt((X1_numerical[0] - X1_analytical[0])**2 + (X1_analytical[1] - X1_numerical[1])**2)


def plotError(lower_N, upper_N):

    n_values = range(lower_N, upper_N)

    errorValuesEuler = [globalError(n, 'euler') for n in n_values]
    errorValuesTrapezoid = [globalError(n, 'trapezoid') for n in n_values]

    plt.loglog(n_values, errorValuesEuler)
    plt.loglog(n_values, errorValuesTrapezoid)
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
    numericalParticleTrajectoryLargeDrag(X0, )
    t_euler = time.time() - t_euler

    t_rk2 = time.time()
    # simpleParticleTrajectory(X0, 0, timeStep_rk2, TIME, rk2)
    t_rk2 = time.time() - t_rk2

    return t_rk2 / t_euler


def main():

    # plotNumericalvsAnalytical(800, 800)
    plotError(800, 1000)
    """

    stepsRequired_euler = halfPointSolver(10, 300, euler, 10)
    stepsRequired_rk2 = halfPointSolver(10, 300, rk2, 10)

    plotNumericalvsAnalytical(stepsRequired_euler, stepsRequired_rk2)

    print('Time rk2/euler = ', timeDifference(stepsRequired_euler, stepsRequired_rk2))
    """
