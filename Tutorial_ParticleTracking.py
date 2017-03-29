import xarray as xr
from Tutorial_ExampleFileFromWiki import Interpolator
import numpy as np


def simulateParticlesMoving():
    dataPath = "NorKyst-800m.nc"
    dataSet = xr.open_dataset(dataPath)

    velocityField = Interpolator(dataSet=dataSet)
    #  test the object
    time = dataSet.time[3]

    X = np.array([-3e6, -1.3e6]).reshape(2, 1)
    print(velocityField)

