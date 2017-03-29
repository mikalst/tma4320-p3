import numpy as np


def euler(X_now, h, time_now, velocityField):
    dt = h / np.timedelta64(1, 's')
    dx_dt = velocityField(X_now, time_now)
    X_next = X_now + dt*dx_dt
    return X_next


# We need to verify that the methods work on known functions
def rk2(X_now, h, time_now, velocityField):
    dt = h / np.timedelta64(1, 's')
    dx_dt = 1/2 * (velocityField(X_now, time_now) + velocityField(X_now, time_now + h))
    X_next = X_now + dt*dx_dt
    return X_next

