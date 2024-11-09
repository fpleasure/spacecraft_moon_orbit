from parametrs import *
from spacecraft_algorithm import spacecraft_algorithm, make_plot

if __name__ == "__main__":
    u_0 = [0, 0, 0, 0, M_0]
    t_0 = 0
    tau = 0.1
    theta_1 = -0.0035
    theta_2 = 0.18
    t_1 = 300
    t_2 = 500
    solution = spacecraft_algorithm(
        u_0, t_0, tau, theta_1, theta_2, t_1, t_2)
    make_plot(solution)
