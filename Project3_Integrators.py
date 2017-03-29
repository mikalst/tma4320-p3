import numpy as np


def euler(X_now, h, time_now, velocityField):
    dt = h / np.timedelta64(1, 's')
    dx_dt = velocityField(X_now, time_now)
    X_next = X_now + dt*dx_dt
    return X_next