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

    # Number of 'particles' used to simulate concentration.
    N_particles = 10000

    # Particles are placed in random positions within given interval
    X0 = np.zeros((2, N_particles))
    X0[0, :] = np.random.uniform(-3.01e6, -2.99e6, N_particles)
    X0[1, :] = np.random.uniform(-1.21e6, -1.19e6, N_particles)

    # Set how time interval
    start_time = 0
    end_time = 24 * 10

    # Calculate particle positinos in end time
    X1 = simulateConcentrationMoving(start_time, end_time, X0)

    # Plot particle positions as a concentration map
    plotConcentrationOnMap(X1)


def simulateConcentrationMoving(hoursAfterDefaultDateStart, hoursAfterDefaultDateStop, initialPositions_X0):

    # Open dataset
    dataSet = xr.open_dataset('NorKyst-800m.nc')

    # Time step is set (Time resolution)
    h = np.timedelta64(3600, 's')

    # Velocity field of sea is calculated
    velocityField = P_Interpolator.Interpolator(dataSet)

    # Current and final time is set
    time_current = np.datetime64('2017-02-01T00:00:00') + np.timedelta64(hoursAfterDefaultDateStart, 'h')
    time_final = time_current + np.timedelta64(hoursAfterDefaultDateStop, 'h')

    # Particle positions in final time is calculated
    current_positions_X1 = P_Uf.particleTrajectory(initialPositions_X0, time_current, h, time_final,
                                                   velocityField, P_Int.rk2, boolReturnOnlyEnd=True)

    return current_positions_X1


def plotConcentrationOnMap(particlePositions_X0):

    # Path of data
    datapath = './NorKyst-800m.nc'

    # Set grid cell size of approximately 800
    bins = max(int(max(particlePositions_X0[0, :]) - min(particlePositions_X0[0, :]))/800,
               int(max(particlePositions_X0[1, :]) - min(particlePositions_X0[1, :]))/800)

    # Calculate 2d histogram of particlepositions with a resolution of ~800m/pixel
    C, x, y = np.histogram2d(particlePositions_X0[0, :], particlePositions_X0[1, :], bins)

    # Make sure we have NetCDF dataset to get projection info
    d = xr.open_dataset(datapath)

    # Set style of plotting
    plt.style.use('bmh')

    # Create a figure object, and an axes instance, with projection info
    fig = plt.figure(figsize=(12, 8))
    ax = plt.axes(projection=ccrs.NorthPolarStereo())

    # Draw land and coastlines
    land_10m = cfeature.NaturalEarthFeature('physical', 'land', '10m', color='#dddddd')
    ax.add_feature(land_10m)
    ax.coastlines(resolution='10m')

    # Convert coordinates to longitude and latitude
    p1 = pyproj.Proj(d.projection_stere.proj4)
    # Next, create latlong projection object
    p2 = pyproj.Proj(proj='latlong')

    lons, lats = pyproj.transform(p1, p2, x, y)

    # Plot data
    ax.pcolormesh(lons, lats, C, transform=ccrs.PlateCarree(), zorder=2)

    # Set the extent of the map.
    ax.set_extent((-5, 15, 57, 67))
    plt.show()
