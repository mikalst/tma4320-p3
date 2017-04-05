"""Main file
"""

import Project3_ParticleTracking as P_Pt
import Project3_ConcentrationTracking as P_Ct
import Project3_TestsWithoutDrag as P_Tc


if __name__ == "__main__":
    MasterFlag = {
        -1: "TestSpace",
        0: "Task1",
        1: "Task2",
        2: "Task3"
    }[1]

    print(MasterFlag)

    if MasterFlag == "TestSpace":
        pass

    elif MasterFlag == "Task1":
        P_Tc.main()

    elif MasterFlag == "Task2":
        P_Pt.main()

    elif MasterFlag == "Task3":
        P_Ct.main()

    else:
        print("Invalid MasterFlag")
