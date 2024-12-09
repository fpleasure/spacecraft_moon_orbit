from parametrs import *
from spacecraft_algorithm import RK4, get_theta, make_plot


if __name__ == "__main__":
    u_0 = [0, 0, 0, 0, M_0]
    t_0 = 0
    tau = 0.1
    theta_1 = -0.004867177409845831
    theta_2 = -0.322283760105
    t_1 = 379.89339053
    t_2 = 766.03681511
    # theta_1, theta_2 = get_theta(
    #     u_0, t_0, tau, theta_1, theta_2, t_1, t_2)
    # print(theta_1, theta_2)
    solution = RK4(u_0, t_0, tau, theta_1, theta_2, t_1, t_2)
    u_n = solution[-1]
    make_plot(solution)
