import numpy as np
from matplotlib import pyplot as plt

m = 1.0E-2  # kg
alpha = 5.0E-5  # Ns/m
L = 1.0e2


def euler(X, t, h):
    koef = 2 * np.pi / (24 * 3600)
    f = np.array([-koef * X[1], koef * X[0]])
    return X + h * f


def global_error(h):
    t_values = np.linspace(0, 48 * 3600, 48 * 3600 / h)
    X = np.array([100, 0])
    for t in t_values:
        X = euler(X, t, h)
    return np.sqrt((X[0] - 100) ** 2 + X[1] ** 2)


def g(h):
    return global_error(h) - 10


def plotGlobalError():
    h_values = [h / 10 for h in range(10, 0, -1)]
    global_error_values = [global_error(h) for h in h_values]

    plt.plot(h_values, global_error_values)
    plt.show()


def halvering(lower_h, upper_h):  # kanskje analytisk?
    TOL = 1E-6

    mid = (lower_h + upper_h) / 2

    if g(mid) > 0:
        upper_h = mid
    else:
        lower_h = mid
    if upper_h - lower_h < TOL:
        return (upper_h + lower_h) / 2
    else:
        print("Foreloepig: ", mid, "\n")
        return halvering(lower_h, upper_h)


def trapezoid(X, t, h):
    euler_step = euler(X, t, h)
    koef = 2 * np.pi / (24 * 3600)
    f = np.array([-koef * X[1], koef * X[0]])
    f_euler_step = np.array([-koef * euler_step[1], koef * euler_step[0]])
    middle = (f + f_euler_step) / 2
    return X + h * middle


def v_w(X):  # water velocity at X
    koef = 2 * np.pi / (24 * 3600)
    return np.array([-koef * X[1], koef * X[0]])


def euler_c(X, v, h):  # v = velocity, a = acceleration
    X += h * v
    a = (alpha / m) * (v_w(X) - v)
    v += h * a
    return X, v


def trapezoid_c(X, v, h):
    X1, v1 = euler_c(X, v, h)
    v = (v1 + v) / 2
    X += h * v
    return X, v


def global_error_c(h):
    return  # loes analytisk, og sammenlign


def get_X_v(X, v, h):  # returnerer euler og trapezoid step og hastighet
    X_euler, v_euler = euler_c(X, v, h)
    v = (v_euler + v) / 2
    X_trapezoid = X + h * v
    return X_trapezoid, X_euler, v


def embedded_pair(X, v, t, h):
    TOL = 100 * h / (48 * 3600)
    X_euler, X_trapezoid, v = get_X_v(X, v, h)
    print(np.linalg.norm(X_euler - X_trapezoid))
    while np.linalg.norm(X_euler - X_trapezoid) > TOL:
        h /= 2
        X_euler, X_trapezoid, v = get_X_v(X, v, h)
    return X_trapezoid, v, t + h


def numericalParticleTrajectoryWithDrag():
    h = 3600  # s
    t = 0
    t_max = 48*3600

    X = np.array([100.0, 0.0])
    v = np.array([0.0, 0.0])

    x_values = [X[0]]
    y_values = [X[1]]

    while t < t_max:
        X, v, t = embedded_pair(X, v, t, h)
        x_values.append(X[0])
        y_values.append(X[1])

    return x_values, y_values


def analyticalParticleTrajectoryWithDrag():
    t_values = np.linspace(0, 48*3600, 48, endpoint=True)

    k = alpha / m
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

    return X[:, 0].real, X[:, 1].real


def plotNumericalAnalyticalWithDrag():
    x_values_num, y_values_num = numericalParticleTrajectoryWithDrag()
    x_values_ana, y_values_ana = analyticalParticleTrajectoryWithDrag()

    plt.plot(x_values_num, y_values_num)
    plt.plot(x_values_ana, y_values_ana)
    plt.legend(['Numerical', 'Analytical'])
    plt.grid(True)
    plt.show()


def main():
    plotNumericalAnalyticalWithDrag()


if __name__ == "__main__":
    main()



