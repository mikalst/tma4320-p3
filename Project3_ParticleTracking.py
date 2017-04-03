import Project3_Integrators as P_Int
import Project3_Interpolator as P_Interpolator  # Long name to differentiate from integrator
import Project3_UtilityFunctions as P_Uf

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import pyproj
import cartopy.crs as ccrs
import cartopy.feature as cfeature


def main():
    simulateParticlesMoving()


def simulateParticlesMoving():

    hoursAfterDefaultDate = range(0, 24*10, 24)
    X0 = np.array([-3e6, -1.2e6]).reshape(2, 1)

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

    plotSimulatedParticlePaths(arrayOf_X1)


def plotSimulatedParticlePaths(array_of_several_X1):
    # plt.style.use('bmh')

    dataSet = xr.open_dataset('NorKyst-800m.nc')

    fig = plt.figure(figsize=(12, 8))
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
    for X1 in array_of_several_X1:
        lons, lats = pyproj.transform(p1, p2, X1[:, 0], X1[:, 1])
        ax.plot(lons, lats, transform=ccrs.PlateCarree(), zorder=2)

    # Set the extent of the map.
    ax.set_extent((-5, 15, 57, 67))

    legend = ['01/02/17', '02/02/17', '03/02/17',
              '04/02/17', '05/02/17', '06/02/17',
              '07/02/17', '08/02/17', '09/02/17',
              '10/02/17']

    plt.legend(legend)

    plt.show()


