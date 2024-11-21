import numpy as np
import matplotlib.pyplot as plt
from typing import Callable

from parametrs import *
from runge_kutta import runge_kutta_4_step, runge_kutta_4_condition
from files import *


def get_secant_method_step_for_theta(u: list, t: float, tau: float, theta_1: float,
                                     theta_2: float, t_1: float, t_2: float,
                                     delta: float = 1e-12):
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
    def function_right_side_1(u, T): return FUNCTION_RIGHT_SIDE(
        u, T, theta_1, theta_2, t_1, t_2)
    u_1 = runge_kutta_4_condition(u, t, tau, function_right_side_1)[-1]

    def function_right_side_2(u, T): return FUNCTION_RIGHT_SIDE(
        u, T, theta_1 + delta, theta_2, t_1, t_2)
    u_2 = runge_kutta_4_condition(u, t, tau, function_right_side_2)[-1]

    def function_right_side_3(u, T): return FUNCTION_RIGHT_SIDE(
        u, T, theta_1, theta_2 + delta, t_1, t_2)
    u_3 = runge_kutta_4_condition(u, t, tau, function_right_side_3)[-1]

    # Вычисляем значения функций условий для вычисления Якобиана
    dtheta_dt_1 = ANGLE_CONDITION(u_1[0], u_1[1], u_1[2], u_1[3])
    dtheta_dt_2 = ANGLE_CONDITION(u_2[0], u_2[1], u_2[2], u_2[3])
    dtheta_dt_3 = ANGLE_CONDITION(u_3[0], u_3[1], u_3[2], u_3[3])

    dr_dt_1 = RADIUS_VECTOR_CONDITION(u_1[2], u_1[3])
    dr_dt_2 = RADIUS_VECTOR_CONDITION(u_2[2], u_2[3])
    dr_dt_3 = RADIUS_VECTOR_CONDITION(u_3[2], u_3[3])

    J = [[(dtheta_dt_2 - dtheta_dt_1) / delta, (dtheta_dt_3 - dtheta_dt_1) / delta],
         [(dr_dt_2 - dr_dt_1) / delta, (dr_dt_3 - dr_dt_1) / delta]]

    print(f"[INFO] Calculated jacobian: J[0, 0] = {
        J[0][0]:.3f}, J[0, 1] = {J[0][1]:.3f}, J[1, 0] = {
        J[1][0]:.3f}, J[1, 1] = {J[1][1]:.3f}.")
    print(f"[INFO] det|J| = {J[0][0] * J[1][1] - J[0][1] * J[1][0]}")

    # Шаг метода секущих
    F = [dtheta_dt_1, dr_dt_1]
    J_inv = []
    try:
        J_inv = np.linalg.inv(J)
    except np.linalg.LinAlgError:
        print("[ERROR] det|J| = 0.")
        return theta_1, theta_2
    delta_J = -J_inv @ F

    # Новые значения theta
    theta_1, theta_2 = theta_1 + delta_J[0], theta_2 + delta_J[1]

    return theta_1, theta_2


def get_theta(u: list, t: float, tau: float, theta_1: float,
              theta_2: float, t_1: float, t_2: float, delta: float = 1e-5,
              eps_radius_vector: float = 1e-7, eps_angle: float = 1e-9):
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
    print("[INFO] Start get_theta.")
    stop_condition_radius_vector = RADIUS_VECTOR_CONDITION(u[2], u[3])
    stop_condition_angle = ANGLE_CONDITION(u[0], u[1], u[2], u[3])

    # Цикл метода секущих
    while abs(stop_condition_radius_vector) > eps_radius_vector or abs(stop_condition_angle) > eps_angle:
        # Делаем шаг метода секущих
        theta_1, theta_2 = get_secant_method_step_for_theta(
            u, t, tau, theta_1, theta_2, t_1, t_2, delta)

        # Считаем с новыми траекторию theta
        def function_right_side(u, T): return FUNCTION_RIGHT_SIDE(
            u, T, theta_1, theta_2, t_1, t_2)
        u_new = runge_kutta_4_condition(u, t, tau, function_right_side)[-1]

        # Обновление значений критерия останова
        stop_condition_radius_vector = RADIUS_VECTOR_CONDITION(
            u_new[2], u_new[3])
        stop_condition_angle = ANGLE_CONDITION(
            u_new[0], u_new[1], u_new[2], u_new[3])
    print("[INFO] Thetas finded.")
    return theta_1, theta_2


def spacecraft_algorithm(u_0: list, t_0: float, tau: float, theta_1: float, theta_2:
                         float, t_1: float, t_2: float, max_iterations: int = 1500000,
                         eps: float = 1e-9, delta: float = 1e-12):
    """Определение программы выведения КА на орбиту искусственного спутника Луны

    Параметры:
    ----------
    - u_0 (list): Начальное состояние космического аппарата в виде списка [Vx, Vy, x, y, m],
      где Vx и Vy — скорости, x и y — координаты, m — масса.
    - t_0 (float): Начальное время.
    - tau (float): Начальный шаг интегрирования.
    - theta_1 (float): theta 1 с точкой.
    - theta_2 (float): theta 2.
    - t_1 (float): Время начала первой фазы управления тягой.
    - t_2 (float): Время начала второй фазы управления тягой.
    - max_iterations (int): Максимальное количество итераций, по умолчанию 15000.
    - eps (float): Точность условия остановки по скорости, по умолчанию 1e-5.
    - delta (float): Точность вычисления производных.

    Возвращает:
    ----------
    - solution (list of lists): Список состояний космического аппарата на каждом временном шаге.
    """
    print("[INFO] Run spacecraft_algorithm.")
    initialize_output_file()
    solution = [u_0[:], ]
    u = u_0[:]
    t = t_0
    iterations = 0
    find_theta = False  # Индикатор нахождения theta_1 и theta_2
    stop_condition = SPEED_CONDITION(u[0], u[1])

    def function_right_side(u, t): return FUNCTION_RIGHT_SIDE(
        u, t, theta_1, theta_2, t_1, t_2)

    while abs(stop_condition) > eps:
        log_iteration_values(t, u)
        print(f"[INFO] Make step spacecraft_algorithm: t: {t}, V_x: {u[0]:.3f}, V_y: {
            u[1]:.3f}, x: {u[2]:.3f}, y: {u[3]:.3f}, m: {u[4]:.3f}.")
        print(f"[INFO] Stop condition (speed condition): {
            stop_condition:.7f}.")
        print(
            f"[INFO] Stop condition (radius-vector): {RADIUS_VECTOR_CONDITION(u[2], u[3]):.15f}")
        print(f"[INFO] Stop condition (angle): {
            ANGLE_CONDITION(u[0], u[1], u[2], u[3]):.7f}")

        # Отыскание theta_1 и theta_2 с помошью метода секущих
        '''
        if ((t > T_V and t <= t_1) or t >= t_2) and not find_theta:
            theta_1, theta_2 = get_theta(
                u, t, tau, theta_1, theta_2, t_1, t_2, delta)
            find_theta = True

            def function_right_side(u, t): return FUNCTION_RIGHT_SIDE(
                u, t, theta_1, theta_2, t_1, t_2)
        '''

        if t < t_1 and t > T_V:
            u = runge_kutta_4_step(u, t, tau, function_right_side)
            if t + tau >= t_1 and (t + tau - t_1) > eps:
                u = runge_kutta_4_step(u, t, t_1 - t, function_right_side)
                log_iteration_values(t_1, u)
            t += tau
            iterations += 1
            continue

        # Шаг метода Рунге-Кутты
        u = runge_kutta_4_step(u, t, tau, function_right_side)

        # Дробление шага
        if stop_condition < eps and SPEED_CONDITION(u[0], u[1]) > eps:
            tau = tau / 2.
            u = solution[-1]
            print(f"[INFO] Crash step: tau: {tau}.")

        stop_condition = SPEED_CONDITION(u[0], u[1])

        if (t < T_V and t + tau > T_V and (t + tau - T_V) > eps):
            u = runge_kutta_4_step(u, t, T_V - t, function_right_side)
            log_iteration_values(T_V, u)

        if (t < t_2 and t + tau > t_2 and (t + tau - t_2) > eps):
            u = runge_kutta_4_step(u, t, t_2 - t, function_right_side)
            log_iteration_values(t_2, u)

        solution.append(u[:])
        t += tau
        iterations += 1

        # Критерий останова по количеству итераций
        if iterations >= max_iterations:
            print("[ERROR] Max iterations!")
            break

    print(f"[INFO] spacecraft_algorithm stopped!")
    return solution


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
    print(f"[INFO] Save plot at '{filename}'.")
