############################################
#  Example code for Project 3 TMA4320 #
############################################
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


###########################################
# Convenience class for reading datafiles #
###########################################

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


###########################################
# Reading velocity data from NetCDF files #
###########################################

# To read the velocity data from the downloaded NetCDF file,
# use the class defined above in the following way:

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
# First, create projection with metadata from dataset
p1 = pyproj.Proj(d.projection_stere.proj4)
# Next, create latlong projection object
p2 = pyproj.Proj(proj='latlong')
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
