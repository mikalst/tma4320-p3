import Project3_Integrators as P_Int
import Project3_Interpolator as P_Interpolator  # Long name to differentiate from integrator
import Project3_UtilityFunctions as P_Uf

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import pyproj
import cartopy.crs as ccrs
import cartopy.feature as cfeature

from decimal import Decimal


def main():

    X0 = np.array([-3.0e6, -1.2e6]).reshape(2, 1)

    arrayOf_X1 = simulateParticlesMoving(X0)
    plotParticlePath(arrayOf_X1, X0)


def simulateParticlesMoving(X0):

    hoursAfterDefaultDate = [24*n for n in range(0, 10, 4)]

    # Open dataset
    dataSet = xr.open_dataset('NorKyst-800m.nc')

    # The continuous velocity field is calculated through spine interpolation, see Project3_Interpolator.
    velocityField = P_Interpolator.Interpolator(dataSet)

    # Set time increment
    h = np.timedelta64(3600, 's')

    # Prepare empty array
    arrayOf_X1 = []

    for start_time in hoursAfterDefaultDate:
        print(start_time)
        time_initial = np.datetime64('2017-02-01T12:00:00') + np.timedelta64(start_time, 'h')
        time_final = time_initial + np.timedelta64(10, 'D')

        # Position vector is calculated from initial position and time
        X1 = P_Uf.particleTrajectory(X0, time_initial, h, time_final, velocityField, P_Int.rk2)
        arrayOf_X1.append(X1)

    return arrayOf_X1


def plotParticlePath(arrayOf_X1, X0):

    # Plot x- and y-values
    plt.figure(0)

    for X1 in arrayOf_X1:
        plt.plot(X1[:, 0], X1[:, 1])

    plt.title(r'$x_0 = {:.2E}$'.format(Decimal(X0[0, 0])) + ', ' + r'$y_0 = {:.2E}$'.format(Decimal(X0[1, 0])))
    plt.xlabel(r'$x$')
    plt.ylabel(r'$y$')
    plt.legend(['01/02/17', '05/02/17', '09/02/17'])
    plt.grid(True)
    plt.savefig('test2a{}.png'.format(str(X0[0][0])), dpi=300)

    # Plot map values
    plt.figure(1)

    # Open map dataset
    dataSet = xr.open_dataset('NorKyst-800m.nc')

    fig = plt.figure(1, figsize=(12, 8))
    ax = plt.axes(projection=ccrs.NorthPolarStereo())

    # Draw land and coastlines, we use built-in functions
    land_10m = cfeature.NaturalEarthFeature('physical', 'land', '10m', color='#dddddd')
    ax.add_feature(land_10m)
    ax.coastlines(resolution='10m')

    # Create projection with metadata from dataset
    p1 = pyproj.Proj(dataSet.projection_stere.proj4)
    # Next, create latlong projection object
    p2 = pyproj.Proj(proj='latlong')
    # Convert coordinates

    # Draw trajectories on the map.
    for X1 in arrayOf_X1:
        lons, lats = pyproj.transform(p1, p2, X1[:, 0], X1[:, 1])
        ax.plot(lons, lats, transform=ccrs.PlateCarree(), zorder=2)

    # Set the extent of the map.
    center_lon, center_lat = pyproj.transform(p1, p2, X0[0], X0[1])
    ax.set_extent((center_lon - 3, center_lon + 4, center_lat - 2, center_lat + 2))

    plt.legend(['01/02/17', '05/02/17', '09/02/17'])
    plt.savefig('test2b{}.png'.format(str(X0[0][0])), dpi=400)
    plt.show()
