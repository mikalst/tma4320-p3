import numpy as np

# 2D spline interpolation routine
from scipy.interpolate import RectBivariateSpline


class Interpolator():
    def __init__(self, dataset):
        self.dataset = dataset

    def get_interpolators(self, X, it):
        # Add a buffer of cells around the extent of the particle cloud
        buf = 3
        # Find extent of particle cloud in terms of indices
        imax = np.searchsorted(self.dataset.X, np.amax(X[0, :])) + buf
        imin = np.searchsorted(self.dataset.X, np.amin(X[0, :])) - buf
        jmax = np.searchsorted(self.dataset.Y, np.amax(X[1, :])) + buf
        jmin = np.searchsorted(self.dataset.Y, np.amin(X[1, :])) - buf
        # Take out subset of array, to pass to
        # interpolation object
        # Fill NaN values (land cells) with 0, otherwise
        # interpolation won't work
        u = self.dataset.u[it, 0, jmin:jmax, imin:imax].fillna(0.0)
        v = self.dataset.v[it, 0, jmin:jmax, imin:imax].fillna(0.0)
        # RectBivariateSpline essentially returns a function,
        # which can be called to get value at arbitrary position
        # kx and ky sets order of spline interpolation along either direction (must be 1 <= kx <= 5)
        # transpose arrays to switch order of coordinates
        fu = RectBivariateSpline(self.dataset.X[imin:imax], self.dataset.Y[jmin:jmax], u.T)  # kx = 3, ky = 3)
        fv = RectBivariateSpline(self.dataset.X[imin:imax], self.dataset.Y[jmin:jmax], v.T)  # kx = 3, ky = 3)
        return fu, fv

    def get_time_index(self, t):
        # Get index of largest timestamp smaller than (or equal to) t
        return np.searchsorted(self.dataset.time, t, side='right') - 1

    def __call__(self, X, t):
        # get index of current time in dataset
        it = self.get_time_index(t)
        # get interpolating functions,
        # covering the extent of the particle
        fu, fv = self.get_interpolators(X, it)
        # Evaluate velocity at position(x[:], y[:])
        dx = fu(X[0, :], X[1, :], grid=False)
        dy = fv(X[0, :], X[1, :], grid=False)
        return np.array([dx, dy])
