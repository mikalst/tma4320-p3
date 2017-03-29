"""Main file
"""

import Project3_ParticleTracking as Tut_Pt

import numpy as np


def oppgave2():

    # Set initial position
    X0 = np.array([-3.1e6, -1.2e6]).reshape(2, 1)

    # Prepare array of position vectors
    array_of_several_X1 = []

    # Find position vectors for all start times
    for start_time in [24*day for day in range(10)]:
        array_of_several_X1.append(Tut_Pt.simulateParticlesMoving(start_time, X0))

    legend = ['01/02/17', '02/02,17', '03/02/17',
              '04/02/17', '05/02,17', '06/02/17',
              '07/02/17', '08/02,17', '09/02/17',
              '10/02/17']

    Tut_Pt.plotSimulatedParticlePaths(array_of_several_X1, legend)


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
        oppgave2()
    else:
        print("Invalid MasterFlag")



