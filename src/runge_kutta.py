import numpy as np
from typing import Callable

from parametrs import *


def runge_kutta_4_step(u: list, t: float, tau: float,
                       function_right_side: Callable):
    """
    Выполняет один шаг метода Рунге-Кутты 4-го порядка для системы ОДУ.

    Параметры:
    ----------
    - u (list): Текущее значение векторного состояния системы.
    - t (float): Текущее время.
    - tau (float): Шаг интегрирования.
    - function_right_side (Callable): Функция правой части ОДУ, которая
      принимает параметры (u, t) и возвращает значение производной в данной точке.

    Возвращает:
    ----------
    - u_next (list): Значение состояния системы на следующем временном шаге.
    """
    k1 = function_right_side(u, t)
    k2 = function_right_side(
        [ui + tau * k1_i / 2 for ui, k1_i in zip(u, k1)], t + tau / 2)
    k3 = function_right_side(
        [ui + tau * k2_i / 2 for ui, k2_i in zip(u, k2)], t + tau / 2)
    k4 = function_right_side(
        [ui + tau * k3_i for ui, k3_i in zip(u, k3)], t + tau)

    u_next = [ui + (tau / 6) * (k1_i + 2 * k2_i + 2 * k3_i + k4_i)
              for ui, k1_i, k2_i, k3_i, k4_i in zip(u, k1, k2, k3, k4)]

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
