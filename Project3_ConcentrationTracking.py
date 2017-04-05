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

    hoursAfterDefaultDate = 24*10

    X1 = calculateParticlePositions(hoursAfterDefaultDate)
    plotParticlePositions(X1, hoursAfterDefaultDate)
    plotConcentration(X1, hoursAfterDefaultDate)


def calculateParticlePositions(end_time):

    # Number of 'particles' used to simulate concentration.
    N_particles = 10000

    # Particles are placed in random positions within given interval
    X0 = np.zeros((2, N_particles))
    X0[0, :] = np.random.uniform(-3.01e6, -2.99e6, N_particles)
    X0[1, :] = np.random.uniform(-1.21e6, -1.19e6, N_particles)

    # Set how time interval
    start_time = 0

    # Open dataset
    dataSet = xr.open_dataset('NorKyst-800m.nc')

    # Time step is set (Time resolution)
    h = np.timedelta64(3600, 's')

    # Velocity field of sea is calculated
    velocityField = P_Interpolator.Interpolator(dataSet)

    # Current and final time is set
    time_current = np.datetime64('2017-02-01T12'
                                 ':00:00') + np.timedelta64(start_time, 'h')
    time_final = time_current + np.timedelta64(end_time, 'h')

    # Particle positions in final time is calculated
    X1 = P_Uf.particleTrajectory(X0, time_current, h, time_final, velocityField, P_Int.rk2, boolReturnOnlyEnd=True)

    return X1


def plotConcentration(X1, hoursAfterDefaultDate):

    # Open dataset
    dataSet = xr.open_dataset('NorKyst-800m.nc')

    # Set style of plotting
    plt.style.use('bmh')

    # Prepare concentration figure
    fig = plt.figure(0, figsize=(12, 8))
    ax = plt.axes(projection=ccrs.NorthPolarStereo())

    # Draw land and coastlines
    land_10m = cfeature.NaturalEarthFeature('physical', 'land', '10m', color='#dddddd')
    ax.add_feature(land_10m)
    ax.coastlines(resolution='10m')

    # Convert coordinates to longitude and latitude
    p1 = pyproj.Proj(dataSet.projection_stere.proj4)
    # Next, create latlong projection object
    p2 = pyproj.Proj(proj='latlong')

    lons, lats = pyproj.transform(p1, p2, X1[0, :], X1[1, :])

    bins = max(int((max(X1[0, :]) - min(X1[0, :]))/800), int((max(X1[1, :]) - min(X1[1, :]))/800))

    # Calculate 2d histogram of particlepositions with a resolution of ~800m/pixel
    C, xedges, yedges = np.histogram2d(lons, lats, bins=bins)

    # Plot data
    im = ax.pcolormesh(xedges, yedges, C.T, transform=ccrs.PlateCarree(), zorder=2, vmax=45) # xedges, yedges

    # Plot colorbar for showing densite unit
    cbar = plt.colorbar(im, ax=ax)
    cbar.ax.set_ylabel('Partikler i gridpunkt, ' + r'$[N/64000m^2]$')

    # Save figure
    time = str(np.datetime64('2017-02-01T12:00:00') + np.timedelta64(hoursAfterDefaultDate, 'h'))
    plt.title(time)
    plt.savefig('concentrationImages/conc{}.png'.format(time))
    plt.show()
    plt.clf()


def plotParticlePositions(X1, hoursAfterDefaultDate):

    # Open dataset
    dataSet = xr.open_dataset('NorKyst-800m.nc')

    # Set style of plotting
    plt.style.use('bmh')

    # Prepare partcile position figure
    fig = plt.figure(1, figsize=(12, 8))
    ax = plt.axes(projection=ccrs.NorthPolarStereo())

    # Draw land and coastlines
    land_10m = cfeature.NaturalEarthFeature('physical', 'land', '10m', color='#dddddd')
    ax.add_feature(land_10m)
    ax.coastlines(resolution='10m')

    # Convert coordinates to longitude and latitude
    p1 = pyproj.Proj(dataSet.projection_stere.proj4)
    p2 = pyproj.Proj(proj='latlong')

    # Transform x- and y-coordinates to longitudes and latitudes
    lons, lats = pyproj.transform(p1, p2, X1[0, :], X1[1, :])

    # Plot data
    ax.scatter(lons, lats, transform=ccrs.PlateCarree(), zorder=2, marker='.')

    ax.set_extent((-2, 8, 58.5, 61))

    time = str(np.datetime64('2017-02-01T12:00:00') + np.timedelta64(hoursAfterDefaultDate, 'h'))
    plt.title(time)
    plt.savefig('part{}.png'.format(time))
    plt.show()
    plt.clf()
