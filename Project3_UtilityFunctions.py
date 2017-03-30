import numpy as np


def particleTrajectory(X_initial, time_initial, h, time_final, velocityField, integrator, boolReturnOnlyEnd=False):
    numberOfTimeSteps = int((time_final - time_initial) / h)
    X = np.zeros((numberOfTimeSteps + 1 + 1, *X_initial.shape))
    #  Må her skrive en forklaring på stjerneoperator
    #  Denne henter ut tupler, men hvorfor gjøre vi dette?
    #  Oversett til engelsk

    X[0, :] = X_initial
    time_now = time_initial

    for step in range(numberOfTimeSteps + 1):
        h = min(h, time_final - time_now)

        X[step + 1, :] = integrator(X[step, :], h, time_now, velocityField)
        time_now += h

    # If only endpoint is needed, as in concentration plotting
    if boolReturnOnlyEnd:
        return X[numberOfTimeSteps, :]

    return X
