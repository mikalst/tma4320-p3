
import Project3_PlottingFunctions as Tut_Pf
import Project3_ParticleTracking as Tut_Pt
import Project3_ExampleFileFromWiki as Tut_Effw
import Project3_Integrators as Tut_Int
import Project3_UtilityFunctions as Tut_Uf

import xarray as xr
import numpy as np


if __name__ == "__main__":
    """Docstring"""
    MasterFlag = {
        -1: "TestSpace",
        0: "CheckOutInterpolator",
        1: "PlotOnMap"
    }[-1]

    print(MasterFlag)

    if MasterFlag == "TestSpace":
        dataSet = xr.open_dataset('NorKyst-800m.nc')
        time_initial = dataSet.time[0]
        h = np.timedelta64(3600, 's')
        time_final = time_initial + np.timedelta64(1, 'D')
        velocityField = Tut_Effw.Interpolator(dataSet)

        X0 = np.array([-3e6, -1.3e6]).reshape(2, 1)

        X1 = Tut_Uf.particleTrajectory(X0, time_initial, h, time_final, velocityField, Tut_Int.euler)

    elif MasterFlag == "CheckOutInterpolator":
        Tut_Pt.simulateParticlesMoving()

    elif MasterFlag == "PlotOnMap":
        pass

    else:
        print("Invalid MasterFlag")

