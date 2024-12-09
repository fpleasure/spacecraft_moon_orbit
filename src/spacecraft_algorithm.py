from parametrs import *
from files import *

import numpy as np
import matplotlib.pyplot as plt
from typing import Callable


def runge_kutta_4_step(
    u: np.ndarray, t: float, tau: float, function_right_side: Callable[[np.ndarray, float], np.ndarray]
) -> np.ndarray:
    """
    Выполняет один шаг метода Рунге-Кутты 4-го порядка для системы ОДУ с использованием NumPy.

    Параметры:
    ----------
    - u (np.ndarray): Текущее значение векторного состояния системы.
    - t (float): Текущее время.
    - tau (float): Шаг интегрирования.
    - function_right_side (Callable): Функция правой части ОДУ, которая
      принимает параметры (u, t) и возвращает значение производной в данной точке.

    Возвращает:
    ----------
    - u_next (np.ndarray): Значение состояния системы на следующем временном шаге.
    """
    # Преобразование входных данных в np.ndarray (если еще не массив)
    u = np.asarray(u, dtype=float)

    k1 = np.asarray(function_right_side(u, t), dtype=float)
    k2 = np.asarray(function_right_side(
        u + tau * k1 / 2., t + tau / 2.), dtype=float)
    k3 = np.asarray(function_right_side(
        u + tau * k2 / 2., t + tau / 2.), dtype=float)
    k4 = np.asarray(function_right_side(u + tau * k3, t + tau), dtype=float)

    u_next = u + (tau / 6.) * (k1 + 2. * k2 + 2. * k3 + k4)
    return u_next


def RK4(u_0: list, t_0: float, tau: float, theta_1: float, theta_2:
        float, t_1: float, t_2: float, max_iterations: int = 1500000,
        eps: float = 1e-7, delta: float = 1e-12, log_info=True):
    if log_info:
        initialize_output_file()
    solution = [u_0[:], ]
    u = u_0[:]
    t = t_0
    iterations = 0
    stop_condition = SPEED_CONDITION(u[0], u[1])
    P_flag = True

    def function_right_side(u, t): return FUNCTION_RIGHT_SIDE(
        u, t, theta_1, theta_2, t_1, t_2, P)

    while abs(stop_condition) > eps:
        if P_flag:
            def function_right_side(u, t): return FUNCTION_RIGHT_SIDE(
                u, t, theta_1, theta_2, t_1, t_2, lambda t, t_1, t_2: 10.27 * 10 ** 3)
        else:
            def function_right_side(u, t): return FUNCTION_RIGHT_SIDE(
                u, t, theta_1, theta_2, t_1, t_2, lambda t, t_1, t_2: 0)

        # if log_info:
            # print(f"[INFO] Make step RK4: t: {t}, V_x: {u[0]:.3f}, V_y: {
            #    u[1]:.3f}, x: {u[2]:.3f}, y: {u[3]:.3f}, m: {u[4]:.3f}.")
            # print(f"[INFO] Stop condition (speed condition): {
            #    stop_condition:.15f}.")
            # print(
            #    f"[INFO] Stop condition (radius-vector): {RADIUS_VECTOR_CONDITION(u[2], u[3]):.15f}")
            # print(f"[INFO] Stop condition (angle): {
            #    ANGLE_CONDITION(u[0], u[1], u[2], u[3]):.15f}")

        if t < t_1 and t > T_V:
            if t + tau > t_1 + eps:
                if abs(t - t_1) < eps:
                    t = t_1
                    iterations += 1
                    P_flag = False
                    continue

                t_p = t + tau
                u = runge_kutta_4_step(u, t, t_1 - t, function_right_side)
                t = t_1
                solution.append(u[:])
                log_iteration_values(t, u, THETA(
                    theta_1, theta_2, t_1, t_2, t), ANGLE_CONDITION(u[0], u[1], u[2], u[3]))
                iterations += 1
                P_flag = False

                if abs(t_1 - (t_p)) > eps:
                    def function_right_side(u, t): return FUNCTION_RIGHT_SIDE(
                        u, t, theta_1, theta_2, t_1, t_2, lambda t, t_1, t_2: 0)

                    u = runge_kutta_4_step(
                        u, t, t_p - t_1, function_right_side)
                    t = t_p
                    solution.append(u[:])
                    log_iteration_values(t, u, THETA(
                        theta_1, theta_2, t_1, t_2, t), ANGLE_CONDITION(u[0], u[1], u[2], u[3]))
                    iterations += 1
                continue

            u = runge_kutta_4_step(u, t, tau, function_right_side)
            solution.append(u[:])
            t += tau
            iterations += 1
            log_iteration_values(t, u, THETA(
                theta_1, theta_2, t_1, t_2, t), ANGLE_CONDITION(u[0], u[1], u[2], u[3]))
            continue

        if t < T_V and t + tau > T_V + eps:
            if abs(t - T_V) < eps:
                t = T_V
                iterations += 1
                P_flag = True
                continue
            t_p = t + tau
            u = runge_kutta_4_step(u, t, T_V - t, function_right_side)
            t = T_V
            solution.append(u[:])
            log_iteration_values(t, u, THETA(
                theta_1, theta_2, t_1, t_2, t), ANGLE_CONDITION(u[0], u[1], u[2], u[3]))
            P_flag = True
            iterations += 1
            if abs(T_V - (t_p)) > eps:
                def function_right_side(u, t): return FUNCTION_RIGHT_SIDE(
                    u, t, theta_1, theta_2, t_1, t_2, lambda t, t_1, t_2: 10.27 * 10 ** 3)

                u = runge_kutta_4_step(
                    u, t, t_p - T_V, function_right_side)
                t = t_p
                solution.append(u[:])
                log_iteration_values(t, u, THETA(
                    theta_1, theta_2, t_1, t_2, t), ANGLE_CONDITION(u[0], u[1], u[2], u[3]))
                iterations += 1
            continue

        if t < t_2 and t + tau > t_2:
            if abs(t - t_2) < eps:
                t = t_2
                iterations += 1
                P_flag = True
                continue
            t_p = t + tau
            u = runge_kutta_4_step(u, t, t_2 - t, function_right_side)
            t = t_2
            solution.append(u[:])
            log_iteration_values(t, u, THETA(
                theta_1, theta_2, t_1, t_2, t), ANGLE_CONDITION(u[0], u[1], u[2], u[3]))
            P_flag = True
            iterations += 1

            if abs(t_2 - (t_p)) > eps:
                def function_right_side(u, t): return FUNCTION_RIGHT_SIDE(
                    u, t, theta_1, theta_2, t_1, t_2, lambda t, t_1, t_2: 10.27 * 10 ** 3)
                u = runge_kutta_4_step(
                    u, t, t_p - t_2, function_right_side)
                t = t_p
                solution.append(u[:])
                log_iteration_values(t, u, THETA(
                    theta_1, theta_2, t_1, t_2, t), ANGLE_CONDITION(u[0], u[1], u[2], u[3]))
                iterations += 1
            continue

        # Шаг метода Рунге-Кутты
        u = runge_kutta_4_step(u, t, tau, function_right_side)

        # Дробление шага
        if stop_condition < eps and SPEED_CONDITION(u[0], u[1]) > eps:
            tau = tau / 2.
            t -= tau
            u = solution[-1]
            # if log_info:
            # print(f"[INFO] Crash step: tau: {tau}.")

        stop_condition = SPEED_CONDITION(u[0], u[1])
        solution.append(u[:])
        t += tau
        iterations += 1
        log_iteration_values(t, u, THETA(
            theta_1, theta_2, t_1, t_2, t), ANGLE_CONDITION(u[0], u[1], u[2], u[3]))

        # Критерий останова по количеству итераций
        if iterations >= max_iterations:
            # print("[ERROR] Max iterations!")
            break

    return solution


def get_secant_method_step_for_theta(u: list, t: float, tau: float, theta_1: float,
                                     theta_2: float, t_1: float, t_2: float,
                                     delta: float = 1e-8, P_flag=True):
    """
    Выполняет один шаг метода секущих для определения новых значений параметров theta_1 и theta_2

    Параметры:
    ----------
    - u (list): Вектор текущего состояния системы.
    - t (float): Текущее время.
    - tau (float): Шаг интегрирования по времени.
    - theta_1, theta_2 (float): Текущие значения параметров theta_1, theta_2.
    - t_1, t_2 (float): Параметры времени.
    - delta (float): Точность вычисления производных.

    Возвращает:
    ----------
    - theta_1, theta_2 (float): Обновленные значения параметров theta_1 и theta_2.
    """
    # Считаем РК с приращениями по theta
    u = [0, 0, 0, 0, M_0]
    t = 0
    tau = 0.1
    u_1 = RK4(
        u, t, tau, theta_1, theta_2, t_1, t_2, log_info=False)[-1]

    u = [0, 0, 0, 0, M_0]
    t = 0
    tau = 0.1
    u_2 = RK4(
        u, t, tau, theta_1 + delta, theta_2, t_1, t_2, log_info=False)[-1]

    u = [0, 0, 0, 0, M_0]
    t = 0
    tau = 0.1
    u_3 = RK4(
        u, t, tau, theta_1, theta_2 + delta, t_1, t_2, log_info=False)[-1]

    # Вычисляем значения функций условий для вычисления Якобиана
    dtheta_dt_1 = ANGLE_CONDITION(u_1[0], u_1[1], u_1[2], u_1[3])
    dtheta_dt_2 = ANGLE_CONDITION(u_2[0], u_2[1], u_2[2], u_2[3])
    dtheta_dt_3 = ANGLE_CONDITION(u_3[0], u_3[1], u_3[2], u_3[3])

    dr_dt_1 = RADIUS_VECTOR_CONDITION(u_1[2], u_1[3])
    dr_dt_2 = RADIUS_VECTOR_CONDITION(u_2[2], u_2[3])
    dr_dt_3 = RADIUS_VECTOR_CONDITION(u_3[2], u_3[3])

    J = [[(dtheta_dt_2 - dtheta_dt_1) / delta, (dtheta_dt_3 - dtheta_dt_1) / delta],
         [(dr_dt_2 - dr_dt_1) / delta, (dr_dt_3 - dr_dt_1) / delta]]

    # print(f"[INFO] Calculated jacobian: J[0, 0] = {
    #    J[0][0]: .5f}, J[0, 1] = {J[0][1]: .5f}, J[1, 0] = {
    #    J[1][0]: .5f}, J[1, 1] = {J[1][1]: .5f}.")
    # print(f"[INFO] det|J| = {J[0][0] * J[1][1] - J[0][1] * J[1][0]}")

    # Шаг метода секущих
    F = [dtheta_dt_1, dr_dt_1]
    J_inv = []
    try:
        J_inv = np.linalg.inv(J)
    except np.linalg.LinAlgError:
        # print("[ERROR] det|J| = 0.")
        return theta_1, theta_2
    delta_J = -J_inv @ F

    # Новые значения theta
    theta_1, theta_2 = theta_1 + delta_J[0], theta_2 + delta_J[1]

    return theta_1, theta_2


def get_theta(u: list, t: float, tau: float, theta_1: float,
              theta_2: float, t_1: float, t_2: float, delta: float = 1e-5,
              eps_radius_vector: float = 1e-7, eps_angle: float = 1e-10):
    """
    Определяет параметры theta_1 и theta_2 с помощью метода секущих для обеспечения
    выполнения заданных условий по радиус-вектору и углу.

    Параметры:
    ----------
    - u (list): Начальное состояние системы.
    - t (float): Начальное время.
    - tau (float): Шаг интегрирования по времени.
    - theta_1, theta_2 (float): Начальные значения параметров управления.
    - t_1, t_2 (float): Параметры времени.
    - delta (float): Точность вычисления производных.
    - eps (float): Точность условия остановки по скорости, по умолчанию 1e-5.

    Возвращает:
    ----------
    - theta_1, theta_2 (float): Оптимальные значения параметров theta_1 и theta_2.
    """
    # print("[INFO] Start get_theta.")
    stop_condition_radius_vector = RADIUS_VECTOR_CONDITION(u[2], u[3])
    stop_condition_angle = ANGLE_CONDITION(u[0], u[1], u[2], u[3])
    iterations = 1
    max_iterations = 20

    # Цикл метода секущих
    while abs(stop_condition_radius_vector) > eps_radius_vector or abs(stop_condition_angle) > eps_angle:
        # Делаем шаг метода секущих
        u_0 = [0, 0, 0, 0, M_0]
        t_0 = 0
        tau = 0.1
        theta_1, theta_2 = get_secant_method_step_for_theta(
            u, t, tau, theta_1, theta_2, t_1, t_2, delta)

        u_0 = [0, 0, 0, 0, M_0]
        t_0 = 0
        tau = 0.1
        u_new = RK4(
            u, t, tau, theta_1, theta_2, t_1, t_2, log_info=False)[-1]

        # Обновление значений критерия останова
        stop_condition_radius_vector = RADIUS_VECTOR_CONDITION(
            u_new[2], u_new[3])
        stop_condition_angle = ANGLE_CONDITION(
            u_new[0], u_new[1], u_new[2], u_new[3])

        # print(
        #    f"[INFO] Stop condition (radius-vector): {RADIUS_VECTOR_CONDITION(u_new[2], u_new[3]):.15f}")
        # print(f"[INFO] Stop condition (angle): {
        #    ANGLE_CONDITION(u_new[0], u_new[1], u_new[2], u_new[3]): .15f}")
        # print(f"[INFO] Thetas: {theta_1}, {theta_2}")

        iterations += 1
        if iterations > max_iterations:
            break
    # print("[INFO] Thetas finded.")
    return theta_1, theta_2


def make_plot(solution: list, filename="img/spacecraft_algorithm.pdf"):
    """
    Создаёт и сохраняет график траектории космического аппарата по результатам работы алгоритма.

    Параметры:
    ----------
    - solution (list of lists): Список состояний космического аппарата на каждом временном шаге.
    - filename (str): Имя файла, в который будут записан график.

    Сохраняет:
    ----------
    - Файл "img/spacecraft_algorithm.pdf": График траектории движения космического аппарата.
    """
    x = list()
    y = list()

    for step in solution:
        x.append(step[2])
        y.append(step[3])

    plt.plot(x, y)
    plt.xlabel("x, м")
    plt.ylabel("y, м")
    plt.grid()
    plt.savefig(filename)
    # print(f"[INFO] Save plot at '{filename}'.")
