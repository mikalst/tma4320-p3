import xarray as xr
import numpy as np
from scipy.interpolate import RectBivariateSpline
import matplotlib.pyplot as plt
import pyproj
import cartopy.crs as ccrs
import cartopy.feature as cfeature


def rk2(X_now, h, time_now, velocityField):
    dt = h / np.timedelta64(1, 's')
    k1 = velocityField(X_now, time_now)
    k2 = velocityField(X_now + dt*k1, time_now)
    X_next = X_now + dt*0.5*(k1 + k2)
    return X_next


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


class OceanSimulation:

    def __init__(self, startTime, endTime, N):

        # Time step is set (Time resolution)
        self.h = np.timedelta64(3600, 's')
        self.N_steps = (endTime - startTime)
        self.timeStart = np.datetime64('2017-02-01T12:00:00') + np.timedelta64(startTime, 'h')
        self.timeFinal = self.timeStart + np.timedelta64(endTime, 'h')

        # Number of 'particles' used to simulate concentration.
        self.N_particles = N

        # Particles are placed in random positions within given interval
        self.X0 = np.zeros((2, self.N_particles))
        self.X0[0, :] = np.random.uniform(-3.01e6, -2.99e6, self.N_particles)
        self.X0[1, :] = np.random.uniform(-1.21e6, -1.19e6, self.N_particles)

        self.X = np.zeros((self.N_steps + 1 + 1, *(self.X0).shape))

        # Open dataset
        self.dataSet = xr.open_dataset('NorKyst-800m.nc')

        # Velocity field of sea is calculated
        self.velocityField = Interpolator(self.dataSet)

        # Set numerical integration method
        self.integrator = rk2

        return

    def calculateParticleMovement(self):

        self.X[0, :] = self.X0

        timeNow = self.timeStart

        for step in range(self.N_steps + 1):
            h = min(self.h, self.timeFinal - timeNow)

            self.X[step + 1, :] = self.integrator(self.X[step, :], h, timeNow, self.velocityField)
            timeNow += h

        return

    def plotConcentration(self, time):

        time = np.datetime64('2017-02-01T12:00:00') + np.timedelta64(time, 'h')

        # Set style of plotting
        plt.style.use('bmh')

        # Prepare concentration figure
        fig = plt.figure(figsize=(8, 8))
        ax = plt.axes(projection=ccrs.NorthPolarStereo())

        # Draw land and coastlines
        land_10m = cfeature.NaturalEarthFeature('physical', 'land', '10m', color='#dddddd')
        ax.add_feature(land_10m)
        ax.coastlines(resolution='10m')

        # Convert coordinates to longitude and latitude
        p1 = pyproj.Proj(self.dataSet.projection_stere.proj4)
        # Next, create latlong projection object
        p2 = pyproj.Proj(proj='latlong')

        # Assert that time to be plotted is already calculated
        assert(self.timeStart <= time <= self.timeFinal)
        X1 = self.X[int((time - self.timeStart) / (self.timeFinal - self.timeStart) * self.N_steps)]

        lons, lats = pyproj.transform(p1, p2, X1[0, :], X1[1, :])

        # Calculate 2d histogram of particlepositions with a resolution of ~800m/pixel
        C, xedges, yedges = np.histogram2d(lons, lats, bins=40)

        # Plot data, set scale: vmax=45
        im = ax.pcolormesh(xedges, yedges, C.T, transform=ccrs.PlateCarree(), zorder=2)  # xedges, yedges,

        # Plot colorbar for showing densite unit
        cbar = plt.colorbar(im, ax=ax)
        cbar.ax.set_ylabel('Partikler i gridpunkt, ' + r'$[N/64000m^2]$')

        # Plot decoration
        plt.title(time)
        plt.savefig('concentrationImages/conc{}.png'.format(time))
        plt.show()
        plt.clf()

    def plotParticlePositions(self, time):

        time = np.datetime64('2017-02-01T12:00:00') + np.timedelta64(time, 'h')

        # Set style of plotting
        plt.style.use('bmh')

        # Prepare partcile position figure
        fig = plt.figure(figsize=(8, 8))
        ax = plt.axes(projection=ccrs.NorthPolarStereo())

        # Draw land and coastlines
        land_10m = cfeature.NaturalEarthFeature('physical', 'land', '10m', color='#dddddd')
        ax.add_feature(land_10m)
        ax.coastlines(resolution='10m')

        # Convert coordinates to longitude and latitude
        p1 = pyproj.Proj(self.dataSet.projection_stere.proj4)
        p2 = pyproj.Proj(proj='latlong')

        assert(self.timeStart <= time <= self.timeFinal)
        X1 = self.X[int((time - self.timeStart) / (self.timeFinal - self.timeStart) * self.N_steps)]

        # Transform x- and y-coordinates to longitudes and latitudes
        lons, lats = pyproj.transform(p1, p2, X1[0, :], X1[1, :])

        # Plot data
        ax.scatter(lons, lats, transform=ccrs.PlateCarree(), zorder=2, marker='.')

        ax.set_extent((-2, 8, 58.5, 61))

        # Plot decoration
        plt.title(time)
        plt.savefig('particleImages/part{}.png'.format(time))
        plt.show()
        plt.clf()


if __name__ == "__main__":
    o = OceanSimulation(0, 24*10, 10000)

    o.calculateParticleMovement()
    o.plotConcentration(24*10)
    o.plotParticlePositions(24*10)
