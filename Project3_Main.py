"""Main file
"""


import Project3_ParticleTracking as Tut_Pt

import numpy as np


if __name__ == "__main__":
    """Docstring"""
    MasterFlag = {
        -1: "TestSpace",
        0: "CheckOutInterpolator",
        1: "PlotOnMap"
    }[1]

    if MasterFlag == "TestSpace":
        pass

    elif MasterFlag == "CheckOutInterpolator":
        pass

    elif MasterFlag == "PlotOnMap":
        X0 = np.array([-3.1e6, -1.2e6]).reshape(2, 1)
        Tut_Pt.simulateParticlesMoving(24, X0=X0)

    else:
        print("Invalid MasterFlag")

