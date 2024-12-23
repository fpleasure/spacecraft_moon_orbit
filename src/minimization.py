import numpy as np
from spacecraft_algorithm import RK4, get_theta
from parametrs import *


def numerical_gradient_and_hessian(func, x, epsilon=1e-4):
    x0, x1 = x
    f0 = func([x0, x1])
    f1 = func([x0 + epsilon, x1])
    f2 = func([x0 - epsilon, x1])
    f3 = func([x0, x1 + epsilon])
    f4 = func([x0, x1 - epsilon])
    f5 = func([x0 + epsilon, x1 + epsilon])

    grad = [(f1 - f2) / (2 * epsilon), (f3 - f4) / (2 * epsilon)]
    h1 = (f1 - 2 * f0 + f2) / (epsilon * epsilon)
    h2 = (f3 - 2 * f0 + f4) / (epsilon * epsilon)
    h = (f5 - f1 - f3 + f0) / (epsilon * epsilon)
    hessian = [[h1, h], [h, h2]]
    return grad, hessian


def numerical_gradient(func, x, epsilon=1e-2):
    grad = np.zeros(2)

    x0, x1 = x
    f0 = func([x0, x1])
    f1 = func([x0 + epsilon, x1])
    f2 = func([x0, x1 + epsilon])

    # Параметры для центральной разности
    x0, x1 = x

    # Производная по x0
    grad[0] = (f1 - f0) / epsilon

    # Производная по x1
    grad[1] = (f2 - f0) / epsilon

    return grad


def numerical_hessian(func, x, epsilon=1e-2):
    hessian = np.zeros((2, 2))

    # Параметры для центральной разности
    x0, x1 = x

    ff = func([x0, x1])

    # Диагональные элементы (вторая производная по x0 и x1)
    hessian[0, 0] = (func([x0 + epsilon, x1]) - 2 *
                     ff + func([x0 - epsilon, x1])) / (epsilon ** 2)
    hessian[1, 1] = (func([x0, x1 + epsilon]) - 2 *
                     ff + func([x0, x1 - epsilon])) / (epsilon ** 2)

    # Вне диагональные элементы (смешанные производные)
    f_pp = func([x0 + epsilon, x1 + epsilon])
    f_pm = func([x0 + epsilon, x1 - epsilon])
    f_mp = func([x0 - epsilon, x1 + epsilon])
    f_mm = func([x0 - epsilon, x1 - epsilon])
    hessian[0, 1] = hessian[1, 0] = (
        f_pp - f_pm - f_mp + f_mm) / (4 * epsilon ** 2)

    return hessian


def gradient_descent_adaptive(func, x0, tol=1e-2, max_iter=25, a0=100, C1=0.5, C2=2):
    x = np.array(x0, dtype=float)
    step_size = a0
    history = [x.copy()]
    func_old = func(x)
    func_new = func_old
    for iteration in range(max_iter):
        # Вычисление градиента
        grad = numerical_gradient(func, x)
        grad_norm = np.linalg.norm(grad)

        # Проверяем, достаточно ли мал градиент
        if grad_norm < tol:
            print(f"[INFO] Метод завершён на итерации {
                  iteration}, ||grad|| = {grad_norm:.6f}")
            break

        # Обновление точки
        x_new = x - step_size * grad
        func_old = func_new
        func_new = func(x_new, print_theta=True)

        while func_new > func_old:
            x_new = x + step_size * grad
            step_size *= C1
            x_new = x - step_size * grad
            func_old = func_new
            func_new = func(x_new, print_theta=True)
            print(f"[INFO] Итерация {
                iteration}: неудачный шаг, шаг уменьшен до {step_size:.6f}")

        step_size *= C2  # Увеличиваем шаг
        x = x_new
        print(f"[INFO] Итерация {
            iteration}: удачный шаг, шаг увеличен до {step_size:.6f}")

        history.append(x.copy())
        print(f"[INFO] x = {x}, func(x) = {
              func(x):.6f}, ||grad|| = {grad_norm:.6f}")

    return x, history


def damped_newton_adaptive(func, x0, tol=1e-2, max_iter=100, a0=1.0, C1=0.5, C2=2):
    x = np.array(x0, dtype=float)
    step_size = a0
    history = [x.copy()]
    func_old = func(x)
    func_new = func_old
    for iteration in range(max_iter):
        # Вычисление градиента
        grad, hessian = numerical_gradient_and_hessian(func, x)
        print(f"[INFO] Count gradient {grad}")
        print(f"[INFO] Count hessian {hessian}")
        grad_norm = np.linalg.norm(grad)

        # Проверяем, достаточно ли мал градиент
        if grad_norm < tol:
            print(f"[INFO] Метод завершён на итерации {
                  iteration}, ||grad|| = {grad_norm:.6f}")
            break

        # Обновление точки
        try:
            delta_x = -np.linalg.solve(hessian +
                                       step_size * np.eye(len(x)), grad)
        except np.linalg.LinAlgError:
            print("[ERROR] Матрица Гессе сингулярна. Метод остановлен.")
            break
        x_new = x + delta_x
        func_old = func_new
        func_new = func(x_new, print_theta=True)

        while func_new > func_old:
            x_new = x - delta_x
            step_size *= C1
            try:
                delta_x = - \
                    np.linalg.solve(hessian + step_size * np.eye(len(x)), grad)
            except np.linalg.LinAlgError:
                print("[ERROR] Матрица Гессе сингулярна. Метод остановлен.")
                break
            x_new = x + delta_x
            func_old = func_new
            func_new = func(x_new, print_theta=True)
            print(f"[INFO] Итерация {
                iteration}: неудачный шаг, шаг уменьшен до {step_size:.6f}")

        step_size *= C2  # Увеличиваем шаг
        x = x_new
        print(f"[INFO] Итерация {
            iteration}: удачный шаг, шаг увеличен до {step_size:.6f}")

        history.append(x.copy())
        print(f"[INFO] x = {x}, func(x) = {
              func(x):.6f}, ||grad|| = {grad_norm:.6f}")

    return x, history


def damped_newton_method(func, x0, tol=1e-2, max_iter=40, a0=100, C1=0.5, C2=None):
    if C2 is None:
        C2 = 1 / C1

    x = np.array(x0, dtype=float)
    a_i = a0
    iter_history = [x.copy()]
    g_min = func(x)
    x_opt = x.copy()
    alpha = 10

    for _ in range(max_iter):

        g, H = numerical_gradient_and_hessian(func, x)
        # g = numerical_gradient(func, x)
        print(f"[INFO] Count gradient {g}")

        if np.linalg.norm(g) < tol:
            break

        # H = numerical_hessian(func, x)
        print(f"[INFO] Count hessian {H}")
        H = H + [[alpha, 0], [0, alpha]]
        Hinv = np.linalg.inv(H)
        try:
            # delta_x = -np.linalg.solve(H + a_i * np.eye(len(x)), g)
            delta_x = -np.linalg.solve(H, g)
        except np.linalg.LinAlgError:
            print("[ERROR] Матрица Гессе сингулярна. Метод остановлен.")
            break

        x_new = x + delta_x
        g_new = func(x_new)

        print(f"Текущая точка: x = {x}")
        print(f"Значение функции: func(x) = {g_new}")
        print(f"Значение функции: func_min(x) = {g_min}")
        print(f"Шаговое изменение: delta_x = {delta_x}")
        print(f"Новое значение: x_new = {x_new}")

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
    print(f"[INFO] H = {H}")

    def f(x, print_theta=False):
        u_0 = [0, 0, 0, 0, M_0]
        t_0 = 0
        tau = 0.1
        t_1 = x[0]
        t_2 = x[1]
        theta_1, theta_2 = -0.004, -0.32
        theta_1, theta_2 = get_theta(
            u_0, t_0, tau, theta_1, theta_2, t_1, t_2)
        solution = RK4(u_0, t_0, tau, theta_1, theta_2, t_1, t_2)
        m_n = solution[-1][-1]
        if print_theta:
            print(f"[INFO] Thetas: {theta_1}, {theta_2}")
        return M_0 - m_n

    # def f(x, print_theta=False):
        # return (1 - x[0]) ** 2 + 100 * (x[1] - x[0] ** 2) ** 2
      #   return x[0] ** 2 + x[1] ** 2
    # Начальная точка
    x0 = [350, 750]

    # Запуск метода
    # opt_x, history = damped_newton_method(f, x0, a0=100)
    # opt_x, history = damped_newton_adaptive(f, x0, a0=100)
    opt_x, history = gradient_descent_adaptive(f, x0, a0=100)

    print("Оптимальная точка:", opt_x)
    print("История итераций:", history)
