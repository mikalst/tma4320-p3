"""Docstring"""

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
    # Set initial position
    X0 = np.array([-3.1e6, -1.2e6]).reshape(2, 1)

    # Prepare array of position vectors
    array_of_several_X1 = []

    # Find position vectors for all start times
    for start_time in [24*day for day in range(10)]:
        array_of_several_X1.append(simulateParticlesMoving(start_time, X0))

    legend = ['01/02/17', '02/02/17', '03/02/17',
              '04/02/17', '05/02/17', '06/02/17',
              '07/02/17', '08/02/17', '09/02/17',
              '10/02/17']

    plotSimulatedParticlePaths(array_of_several_X1, legend)


def simulateParticlesMoving(hoursAfterDefaultDate=0, X0=np.array([-3e6, -1.2e6]).reshape(2, 1)):
    """
    Function used in Task 2, used for only 1 particle with start and end specified in position
    """
    # plt.style.use('bmh')

    dataSet = xr.open_dataset('NorKyst-800m.nc')

    # The initial time is set to now
    time_initial = dataSet.time[hoursAfterDefaultDate]
    h = np.timedelta64(3600, 's')

    # Final time is set to initial time plus a time period of 1 'D' = 1 day. 3600 's' = 3600 seconds, 24 'h' = 24 hours
    time_final = time_initial + np.timedelta64(1, 'D')

    # The continuous velocity field is calculated through spine interpolation, see Project3_Interpolator.
    velocityField = P_Interpolator.Interpolator(dataSet)

    # Secure that initial position is given and valid. If not, set to default value
    # X0 = np.array([-3e6, -1.3e6]).reshape(2, 1)

    # Position vector is calculated from initial position
    X1 = P_Uf.particleTrajectory(X0, time_initial, h, time_final, velocityField, P_Int.rk2)

    return X1


def plotSimulatedParticlePaths(array_of_several_X1, legend=[""]):

    dataSet = xr.open_dataset('NorKyst-800m.nc')

    fig = plt.figure(figsize=(12, 8))
    ax = plt.axes(projection=ccrs.NorthPolarStereo())

    # Step 3:
    # It doesn't look like a map unless we add land and sea
    # In order to draw land and coastlines, we use built-in functions
    # in Cartopy.
    land_10m = cfeature.NaturalEarthFeature('physical', 'land', '10m', color='#dddddd')
    ax.add_feature(land_10m)
    ax.coastlines(resolution='10m')

    # Step 4:
    # Coordinate transformation
    # In order to plot trajectories, we need to convert coordinates
    # from xy to longitude and latitude
    # For this we will use the pyproj library
    # First, create projection with metadata from dataset
    p1 = pyproj.Proj(dataSet.projection_stere.proj4)
    # Next, create latlong projection object
    p2 = pyproj.Proj(proj='latlong')
    # Convert coordinates

    # Step 5:
    # Draw trajectories on the map. The transform argument
    # tells Cartopy how to convert from long-lat to projection
    # coordinates. If you leave out zorder, the lines may be
    # hidden behind the map features.

    for X1 in array_of_several_X1:
        lons, lats = pyproj.transform(p1, p2, X1[:, 0], X1[:, 1])
        ax.plot(lons, lats, transform=ccrs.PlateCarree(), zorder=2)

    # Step 6 (optional):
    # Set the extent of the map. If we leave out these, it would
    # just cover the plotted points, and nothing more. Specify
    # (lon0, lon1, lat0, lat1), and Cartopy will make sure the
    # map area is large enough to cover the four points
    # (lon0, lat0), (lon0, lat1), (lon1, lat0), (lon1, lat1).
    ax.set_extent((-5, 15, 57, 67))

    plt.legend(legend)

    plt.show()


