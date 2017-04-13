import numpy as np


def particleTrajectory(X_initial, time_initial, h, time_final, velocityField, integrator, boolReturnOnlyEnd=False):
    numberOfTimeSteps = int((time_final - time_initial) / h)

    # Prepare an empty array to accomodate all position tuples. The * notation fetches out the tuple of the original
    #   position vector shape. This makes the function general for N particles.
    X = np.zeros((numberOfTimeSteps + 1 + 1, *X_initial.shape))

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
