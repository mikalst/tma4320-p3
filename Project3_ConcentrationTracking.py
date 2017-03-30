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

    N_particles = 50

    X0_x = np.random.uniform(-3.01e6, -2.99e6, N_particles)
    X0_y = np.random.uniform(-1.21e6, -1.19e6, N_particles)

    X0 = np.column_stack((X0_x, X0_y))
    X0 = np.array(X0).reshape(N_particles, 2, 1)

    X0 = np.zeros((2, N_particles))
    X0[0, :] = np.random.uniform(-3.01e6, -2.99e6, N_particles)
    X0[1, :] = np.random.uniform(-1.21e6, -1.19e6, N_particles)

    # Prepare array of position vectors
    start_time = 0
    end_time = 24 * 10

    X1 = simulateConcentrationMoving(start_time, end_time, X0)

    plotConcentrationOnMap(X1)


def simulateConcentrationMoving(hoursAfterDefaultDateStart, hoursAfterDefaultDateStop, initialPositions_X0):
    """
    Function used in task 3, used for N particles with start, end and step specified in arguments. Function yields a
    N*2*1 position array for each step.
    """
    # plt.style.use('bmh')

    dataSet = xr.open_dataset('NorKyst-800m.nc')

    # The initial and current is set to February 1st

    h = np.timedelta64(12*3600, 's')

    velocityField = P_Interpolator.Interpolator(dataSet)

    time_current = np.datetime64('2017-02-01T12:00:00')
    time_final = time_current + np.timedelta64(hoursAfterDefaultDateStop, 'h')
    current_positions_X1 = P_Uf.particleTrajectory(initialPositions_X0, time_current, h, time_final,
                                                   velocityField, P_Int.rk2, boolReturnOnlyEnd=True)

    return current_positions_X1


def plotConcentrationOnMap(particlePositions_X0, positionCenter=[-3e6, -1.2e6]):

    # Path of data
    datapath = './NorKyst-800m.nc'

    # the x and y coordinates, and values of your grid
    Nx, Ny = 200, 200

    cellWidth = 800

    # Define x and y coordinates of grid cells
    x = positionCenter[0] + cellWidth * np.arange(Nx)
    y = positionCenter[1] + cellWidth * np.arange(Ny)

    # Convert grid coordinates to 2d arrays
    x, y = np.meshgrid(x, y)
    # Define grid values
    # Er det her mye å tjene på å bruke np.arrays istedenfor?
    C = x*0 + y*0

    # If particle detected in grid, add concentration
    for X0 in particlePositions_X0:
        pos_x = int((X0[0] - positionCenter[0])/cellWidth) + 100
        pos_y = int((X0[1] - positionCenter[1])/cellWidth) + 100
        C[pos_x][pos_y] += 1

    # Assure that program breaks if not all particles are placed in grid
    assert(int(sum(sum(C))) == len(particlePositions_X0))

    # Step 1:
    # Make sure we have NetCDF dataset to get projection info
    # Adjust path to where you saved file
    d = xr.open_dataset(datapath)

    # Step 2:
    # Create a figure object, and an axes instance, with projection info
    fig = plt.figure(figsize=(12, 8))
    ax = plt.axes(projection=ccrs.NorthPolarStereo())

    # Step 3:
    # It doesn't look like a map unless we add land and sea
    # In order to draw land and coastlines, we use built-in functions
    # in Cartopy. These will download some data from
    # www.naturalearthdata.com/ the first time they are called.
    # (resolution 10m means 1 : 10.000.000, not 10 meters)
    land_10m = cfeature.NaturalEarthFeature('physical', 'land', '10m', color='#dddddd')
    ax.add_feature(land_10m)
    ax.coastlines(resolution='10m')

    # Step 4:
    # Convert coordinates to longitude and latitude
    # First, create projection with metadata from dataset
    p1 = pyproj.Proj(d.projection_stere.proj4)
    # Next, create latlong projection object
    p2 = pyproj.Proj(proj='latlong')

    lons, lats = pyproj.transform(p1, p2, x, y)

    # Step 5:
    # Plot data
    ax.pcolormesh(lons, lats, C, transform=ccrs.PlateCarree(), zorder=2)

    # Step 6 (optional):
    # Set the extent of the map. If we leave out these, it would
    # just cover the plotted points, and nothing more. Specify
    # (lon0, lon1, lat0, lat1), and Cartopy will make sure the
    # map area is large enough to cover the four points
    # (lon0, lat0), (lon0, lat1), (lon1, lat0), (lon1, lat1).
    ax.set_extent((-5, 15, 57, 67))
    plt.show()

