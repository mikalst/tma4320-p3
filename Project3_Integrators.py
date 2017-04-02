import numpy as np


def euler(X_now, h, time_now, velocityField):
    dt = h / np.timedelta64(1, 's')
    dx_dt = velocityField(X_now, time_now)
    X_next = X_now + dt*dx_dt
    return X_next


def rk2(X_now, h, time_now, velocityField):
    dt = h / np.timedelta64(1, 's')
    k1 = velocityField(X_now, time_now)
    k2 = velocityField(X_now + dt*k1, time_now)
    X_next = X_now + dt*0.5*(k1 + k2)
    return X_next
