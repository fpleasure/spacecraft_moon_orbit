import numpy as np
from typing import Callable

from parametrs import *


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


def runge_kutta_4(u_0: list, t_0: float, t_k: float,
                  tau: float, function_right_side: Callable):
    """
    Решает систему ОДУ с помощью метода Рунге-Кутты 4-го порядка.

    Параметры:
    ----------
    - u_0 (list): Начальное значение состояния системы.
    - t_0 (float): Начальное время.
    - t_k (float): Конечное время.
    - tau (float): Шаг интегрирования.
    - function_right_side (Callable): Функция правой части ОДУ, которая
      принимает параметры (u, t) и возвращает значение производной в данной точке.

    Возвращает:
    ----------
    - solution (list of lists): Список состояний системы на каждом временном шаге.
    """
    solution = [u_0[:], ]
    u = u_0[:]
    t = t_0

    while t < t_k:
        u = runge_kutta_4_step(u, t, tau, function_right_side)
        solution.append(u[:])
        t += tau

    return solution


def runge_kutta_4_condition(u_0: list, t_0: float, tau: float,
                            function_right_side: Callable,
                            eps_speed: float = 1e-4,
                            eps_radius_vector: float = 1e-3,
                            eps_angle: float = 1e-5,
                            max_iterations: int = 15000):
    """
    Решает систему ОДУ с помощью метода Рунге-Кутты 4-го порядка.

    Параметры:
    ----------
    - u_0 (list): Начальное значение состояния системы.
    - t_0 (float): Начальное время.
    - tau (float): Шаг интегрирования.
    - function_right_side (Callable): Функция правой части ОДУ, которая
      принимает параметры (u, t) и возвращает значение производной в данной точке.
    - eps (float): очность условия остановки по скорости, по умолчанию 1e-5.
    - max_iterations (int): Максимальное количество итераций, по умолчанию 15000.

    Возвращает:
    ----------
    - solution (list of lists): Список состояний системы на каждом временном шаге.
    """
    solution = [u_0[:], ]
    u = u_0[:]
    t = t_0
    iterations = 0
    stop_condition_speed = SPEED_CONDITION(u[0], u[1])
    stop_condition_radius_vector = RADIUS_VECTOR_CONDITION(u[2], u[3])
    stop_condition_angle = ANGLE_CONDITION(u[0], u[1], u[2], u[3])

    while abs(stop_condition_radius_vector) > eps_radius_vector or abs(stop_condition_angle) > eps_angle:
        u = runge_kutta_4_step(u, t, tau, function_right_side)

        # Дробление шага
        if stop_condition_speed < eps_speed and SPEED_CONDITION(u[0], u[1]) > eps_speed:
            tau = tau / 2.
            u = solution[-1]

        stop_condition_speed = SPEED_CONDITION(u[0], u[1])
        stop_condition_radius_vector = RADIUS_VECTOR_CONDITION(u[2], u[3])
        stop_condition_angle = ANGLE_CONDITION(u[0], u[1], u[2], u[3])

        solution.append(u[:])
        t += tau
        iterations += 1

        # Критерий останова по количеству итераций
        if iterations >= max_iterations:
            break

    return solution
