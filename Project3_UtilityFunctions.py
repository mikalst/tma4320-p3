import numpy as np


def particleTrajectory(X_initial, time_initial, h, time_final, velocityField, integrator):
    numberOfTimeSteps = int((time_final - time_initial) / h)
    X = np.zeros((numberOfTimeSteps + 1 + 1, *X_initial.shape))
    #  Vi vil ha X0 som start til en drøss med partikler
    #  Vi får ut shapen til X0
    #  Stjerne er kommando for å hente ut tupler
    #  Altså vil funksjonen vår være generell for vilkårling antall partikler
    #  Oversett til engelsk

    X[0, :] = X_initial
    time_now = time_initial

    for step in range(numberOfTimeSteps + 1):
        h = min(h, time_final - time_now)
        time_now += h

        X[step + 1, :] = integrator(X[step, :], h, time_now, velocityField)

    return X
