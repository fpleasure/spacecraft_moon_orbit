import numpy as np
from spacecraft_algorithm import RK4, get_theta
from parametrs import *


def numerical_gradient(func, x, epsilon=1e-2):
    grad = np.zeros_like(x)
    for i in range(len(x)):
        x1 = x.copy()
        x2 = x.copy()
        x1[i] += epsilon
        x2[i] -= epsilon
        grad[i] = (func(x1) - func(x2)) / (2 * epsilon)
    return grad


def numerical_hessian(func, x, epsilon=1e-2):
    n = len(x)
    hessian = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            x_ij_plus = x.copy()
            x_ij_minus = x.copy()
            x_i_plus = x.copy()
            x_j_plus = x.copy()

            if i == j:
                x_ij_plus[i] += epsilon
                x_ij_minus[i] -= epsilon
                hessian[i, j] = (func(x_ij_plus) - 2 * func(x) +
                                 func(x_ij_minus)) / (epsilon ** 2)
            else:
                x_ij_plus[i] += epsilon
                x_ij_plus[j] += epsilon
                x_i_plus[i] += epsilon
                x_j_plus[j] += epsilon

                hessian[i, j] = (
                    func(x_ij_plus) - func(x_i_plus) - func(x_j_plus) + func(x)
                ) / (epsilon ** 2)
    return hessian


def damped_newton_method(func, x0, tol=1e-2, max_iter=2, a0=0.01, C1=0.5, C2=None):
    if C2 is None:
        C2 = 1 / C1

    x = np.array(x0, dtype=float)
    a_i = a0
    iter_history = [x.copy()]
    g_min = func(x)
    x_opt = x.copy()

    for _ in range(max_iter):
        g = numerical_gradient(func, x)
        print(f"[INFO] Count gradient {g}")

        if np.linalg.norm(g) < tol:
            break

        H = numerical_hessian(func, x)
        print(f"[INFO] Count gradient {H}")

        try:
            delta_x = -np.linalg.solve(H + a_i * np.eye(len(x)), g)
        except np.linalg.LinAlgError:
            print("[ERROR] Матрица Гессе сингулярна. Метод остановлен.")
            break

        x_new = x + delta_x
        g_new = func(x_new)

        if g_new < g_min:
            x = x_new
            g_min = g_new
            a_i *= C1
            x_opt = x.copy()
        else:
            a_i *= C2

        iter_history.append(x.copy())
        print(f"[INFO] t1, t2 = {x}")

    return x_opt, iter_history


if __name__ == "__main__":

    def f(x):
        u_0 = [0, 0, 0, 0, M_0]
        t_0 = 0
        tau = 0.1
        theta_1 = -0.0048761
        theta_2 = -0.3222
        t_1 = x[0]
        t_2 = x[1]
        theta_1, theta_2 = get_theta(
            u_0, t_0, tau, theta_1, theta_2, t_1, t_2)
        solution = RK4(u_0, t_0, tau, theta_1, theta_2, t_1, t_2)
        m_n = solution[-1][-1]
        return M_0 - m_n

    # Начальная точка
    x0 = [379.89339053, 766.03681511]

    # Запуск метода
    opt_x, history = damped_newton_method(f, x0)

    print("Оптимальная точка:", opt_x)
    print("История итераций:", history)
