import numpy as np
from spacecraft_algorithm import RK4, get_theta
from parametrs import *
from minimization import gradient_descent_adaptive, damped_newton_adaptive, damped_newton_method


def get_dm_pg_dm_const():
    dmass = 20
    x0 = [300, 700]
    u_0 = [0, 0, 0, 0, M_0 + dmass]
    t_0 = 0
    tau = 0.1

    def f(x, print_theta=False):
        theta_1 = -0.00486
        theta_2 = -0.32
        t_1 = 300
        t_2 = 700
        theta_1, theta_2 = get_theta(
            u_0, t_0, tau, theta_1, theta_2, t_1, t_2)
        solution = RK4(u_0, t_0, tau, theta_1, theta_2, t_1, t_2)
        m_n = solution[-1][-1]
        if print_theta:
            print(f"[INFO] Thetas: {theta_1}, {theta_2}")
        return M_0 - m_n + dmass

    print(f"[INFO] 1 optimization")
    opt_x_1, history = gradient_descent_adaptive(f, x0, a0=100)
    solution_1 = RK4(u_0, t_0, tau, opt_x_1[0], opt_x_1[1], x0[0], x0[1])[-1]

    dmass = -dmass

    def f(x, print_theta=False):
        theta_1 = -0.004867177409845837
        theta_2 = -0.3222837601056815
        t_1 = 379.89339053
        t_2 = 766.03681511
        theta_1, theta_2 = get_theta(
            u_0, t_0, tau, theta_1, theta_2, t_1, t_2)
        solution = RK4(u_0, t_0, tau, theta_1, theta_2, t_1, t_2)
        m_n = solution[-1][-1]
        if print_theta:
            print(f"[INFO] Thetas: {theta_1}, {theta_2}")
        return M_0 - m_n + dmass

    print(f"[INFO] 2 optimization")
    opt_x_2, history = gradient_descent_adaptive(f, x0, a0=100)
    solution_2 = RK4(u_0, t_0, tau, opt_x_2[0], opt_x_2[1], x0[0], x0[1])[-1]

    fp = solution_1[-1]
    fm = solution_2[-1]

    print(f"[INFO] fp = {fp}, fm = {fm}")
    print(f"[INFO] dm_пг/dm_const = {(fp - fm) / (2 * abs(dmass))}")


if __name__ == "__main__":
    delta_m_const = 10
    delta_w_eff = 10
    m_const = 615

    def f(dmass):
        x0 = [360, 760]
        theta_1 = -0.004867863608316076
        theta_2 = -0.3184399586591739

        def f_inner(x, print_theta=False):
            theta_1 = -0.004867863608316076
            theta_2 = -0.3184399586591739
            u_0 = [0, 0, 0, 0, M_0 + dmass]
            t_0 = 0
            tau = 0.1
            t_1 = x[0]
            t_2 = x[1]
            solution = RK4(u_0, t_0, tau, theta_1, theta_2, t_1, t_2)
            m_n = solution[-1][-1]
            theta_1, theta_2 = get_theta(
                u_0, t_0, tau, theta_1, theta_2, t_1, t_2, dmass=dmass)
            if print_theta:
                print(f"[INFO] Thetas: {theta_1}, {theta_2}")
            return M_0 + dmass - m_n

        # opt_x, _ = gradient_descent_adaptive(f_inner, x0, a0=100)
        # opt_x, history = gradient_descent_adaptive(f_inner, x0, a0=1000)
        # opt_x, history, thetas = damped_newton_method_with_theta(
        #    [theta_1, theta_2], [t_1, t_2])
        u_0 = [0, 0, 0, 0, M_0 + dmass]
        t_0 = 0
        tau = 0.1

        solution = RK4(u_0, t_0, tau, theta_1, theta_2, x0[0], x0[1])
        m_n = solution[-1][-1]

        return m_n - (m_const + dmass)

    fp = f(delta_m_const)
    fm = f(-delta_m_const)
    print(f"[INFO] fp = {fp}, fm = {fm}")
    print(f"[INFO] dm_пг/dm_const = {(fp - fm) / delta_m_const / 2}")

    def g(dw_eff):
        u_0 = [0, 0, 0, 0, M_0]
        t_0 = 0
        tau = 0.1
        theta_1 = -0.0048761
        theta_2 = -0.3222
        t_1 = 379.89339053
        t_2 = 766.03681511

        # opt_x, _, thetas = damped_newton_method_with_theta(
        #    [theta_1, theta_2], [t_1, t_2])
        # t_1, t_2 = opt_x
        # theta_1, theta_2 = thetas

        solution = RK4(u_0, t_0, tau, theta_1, theta_2,
                       t_1, t_2, W_EFF=3510 + dw_eff)
        m_n = solution[-1][-1]

        return m_n - (m_const)

    gp = g(delta_w_eff)
    gm = g(-delta_w_eff)
    print(f"[INFO] gp = {gp}, gm = {gm}")
    print(f"[INFO] dm/dw = {(gp - gm) / delta_w_eff / 2}")
