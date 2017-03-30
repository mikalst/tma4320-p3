############################################
#  Example code for Project 3 TMA4320 #
############################################

from Project3_Interpolator import Interpolator

import numpy as np
from matplotlib import pyplot as plt
# nicer looking default plots
plt.style.use('bmh')

# Library to read data from NetCDF files
import xarray as xr



# 2D spline interpolation routine
from scipy.interpolate import RectBivariateSpline

# Map plotting library
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# library for coordinate transformations
import pyproj

#Test 2

###########################################
# Convenience class for reading datafiles #
###########################################


def main():
    ###########################################
    # Reading velocity data from NetCDF files #
    ###########################################

    # To read the velocity data from the downloaded NetCDF file,
    # Initialise interpolator with dataset
    # (adjust path to where you saved the file)
    datapath = './NorKyst-800m.nc'
    # Open the datafile as an xarray dataset
    d = xr.open_dataset(datapath)
    # Create the interpolator object
    f = Interpolator(dataset=d)

    # f is now essentially like a function, that you
    # can pass to your integrator. Call it with the
    # time and position to get a velocity:
    t = np.datetime64('2017-02-01T12:00:00')
    # The reshape is because f expects X to have the shape (2, Np)
    # and in this case, Np = 1
    X = np.array([-3000000, -1300000]).reshape(2, 1)
    print(f(X, t))


    #################
    # Date and time #
    #################

    # Note that as the downloaded current data are
    # historical records of the estimated state of the
    # ocean, we need actual time and date references.
    # For this, we use the datetime64 class in numpy:

    # Set initial time
    # Dataset covers 2017-02-01 00:00 to 2017-02-19 23:00
    t0 = np.datetime64('2017-02-01T12:00:00')
    # note that h also has time units, for convenient
    # calculation of t + h
    h = np.timedelta64(3600, 's')

    # Now we can step forward in time using
    i = 2
    t = t0 + i*h

    # And to get number of seconds in a timestep,
    # to multiply with velocity (which is in m/s),
    # we can use
    h_seconds = h / np.timedelta64(1, 's')


    ##################################
    # Plotting trajectories on a map #
    ##################################

    # The example below shows how to plot a trajectory on a map.
    # Note that this requires metadata from the NetCDF dataset

    # Step 0:
    # Create some example data for plotting.
    # To use this code, replace x and y with
    # the x and y positions of your trajectory
    Nt = 1000
    # Create a nice wavy trajectory
    x = -2800000 + 30000*np.linspace(0, 4*np.pi, Nt)
    y = -1200000 + 10000*np.sin(np.linspace(0, 12*np.pi, Nt))

    # Step 1:
    # Make sure we have NetCDF dataset to get projection info
    # Adjust path to where you saved file
    d = xr.open_dataset('NorKyst-800m.nc')

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
    # Coordinate transformation
    # In order to plot trajectories, we need to convert coordinates
    # from xy to longitude and latitude
    # For this we will use the pyproj library

    # Convert coordinates
    lons, lats = pyproj.transform(p1, p2, x, y)

    # Step 5:
    # Draw trajectories on the map. The transform argument
    # tells Cartopy how to convert from long-lat to projection
    # coordinates. If you leave out zorder, the lines may be
    # hidden behind the map features.
    ax.plot(lons, lats, transform=ccrs.PlateCarree(), zorder=2)

    # Step 6 (optional):
    # Set the extent of the map. If we leave out these, it would
    # just cover the plotted points, and nothing more. Specify
    # (lon0, lon1, lat0, lat1), and Cartopy will make sure the
    # map area is large enough to cover the four points
    # (lon0, lat0), (lon0, lat1), (lon1, lat0), (lon1, lat1).
    ax.set_extent((-5, 15, 57, 67))


    ##################################
    # Plotting gridded data on a map #
    ##################################

    # Step 0:
    # Create some example gridded data for plotting.
    # In this case, we will plot a simple parabola.
    # To use this code, replace x, y and C with
    # the x and y coordinates, and values of your grid
    Nx, Ny = 200, 200
    # Define x and y coordinates of grid cells
    x = -3000000 + 500 * np.arange(Nx)
    y = -1200000 + 500 * np.arange(Ny)
    # Convert grid coordinates to 2d arrays
    # giving coordinates of cell i,j as
    # (x[i,j], y[i,j])
    x, y = np.meshgrid(x, y)
    # Define grid values
    C = (x + 2950000)**2 + (y + 1150000)**2

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

if __name__ == "__main__":
    main()
