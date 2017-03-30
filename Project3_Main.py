"""Main file
"""

import Project3_ParticleTracking as P_Pt
import Project3_ConcentrationTracking as P_Ct


if __name__ == "__main__":
    """Docstring"""
    MasterFlag = {
        -1: "TestSpace",
        0: "Task1",
        1: "Task2",
        2: "Task3"
    }[2]

    print(MasterFlag)

    if MasterFlag == "TestSpace":
        pass

    elif MasterFlag == "Task1":
        pass

    elif MasterFlag == "Task2":
        P_Pt.main()

    elif MasterFlag == "Task3":
        P_Ct.main()

    else:
        print("Invalid MasterFlag")



