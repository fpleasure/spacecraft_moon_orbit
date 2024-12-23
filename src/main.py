from parametrs import *
from spacecraft_algorithm import RK4, get_theta, make_plot

if __name__ == "__main__":
    u_0 = [0, 0, 0, 0, M_0]
    t_0 = 0
    tau = 0.1
    # theta_1 = -0.004867177409845837
    # theta_2 = -0.3222837601056815
    # t_1 = 379.89339053
    # t_2 = 766.03681511
    t_1, t_2 = 374.03405, 752.09775
    theta_1, theta_2 = -0.004, -0.3
    theta_1, theta_2 = get_theta(
        u_0, t_0, tau, theta_1, theta_2, t_1, t_2)
    # theta_1, theta_2 = -0.004853903297341381 - 0.30920320677326957
    print(theta_1, theta_2)
    solution = RK4(u_0, t_0, tau, theta_1, theta_2, t_1, t_2)
    u_n = solution[-1]
    print(u_n)
    make_plot(solution)
