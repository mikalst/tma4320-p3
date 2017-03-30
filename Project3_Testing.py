import numpy as np
import Project3_ParticleTracking as P_Pt


def main():
    # Set initial position
    X0_python_array = [-3.1e6, -1.2e6, -3.1e6, -1.2e6, -3.1e6, -1.1e6, -3.1e6, -1.2e6, -3.1e6, -1.2e6, -3.1e6, -1.3e6,
                       -3.1e6, -1.2e6, -3e6, -1.2e6, -3.1e6, -1.2e6, -3.1e6, -1.2e6,]
    X0 = np.array(X0_python_array).reshape(10, 2, 1)

    # Prepare array of position vectors
    start_time = 0
    end_time = 24*10
    step = 24*2

    position_vector_generator = P_Pt.simulateParticlesMoving(start_time, end_time, step, X0)

    print(position_vector_generator.__next__())

    print(position_vector_generator.__next__())

    print(position_vector_generator.__next__())



if __name__ == "__main__":
    main()