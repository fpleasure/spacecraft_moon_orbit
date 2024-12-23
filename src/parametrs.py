import numpy as np
from typing import Callable

# Параметры варианта
M_0: float = 3000  # Начальная масса СВ [кг]
T_V: float = 14  # Время вертикального участка [с]
R_M: float = 1738 * 10 ** 3  # Радиус луны [м]
MU_M: float = 4.903 * 10 ** 12  # Гравитационная постоянная луны [м^3/c^2]
H: float = 177 * 10 ** 3  # Высота круговой орбиты ИСЛ [м]
# W_EFF: float = 3510  # Эффективная скорость истечения топлива [м/c]


def G_X(x: float, y: float):
    """Проекция ускорения от гравитационного
    поля Луны на ось Ох [м/с^2]"""
    return - MU_M * x / pow((x ** 2. + (R_M + y) ** 2.), 1.5)


def G_Y(x: float, y: float):
    """Проекция ускорения от гравитационного
    поля Луны на ось Оу [м/с^2]"""
    return - MU_M * (R_M + y) / pow((x ** 2. + (R_M + y) ** 2.), 1.5)


def P(t: float, t_1: float, t_2: float):
    """Тяга [Н]"""
    # if t <= t_1 or t > t_2:
    return 10.27 * 10 ** 3
    # else:
    #    return 0


# def BETA(t: float, t_1: float, t_2: float):
#     """Секундный расход топлива [кг/с]"""
#     return P(t, t_1, t_2) / W_EFF


def THETA(theta_1: float, theta_2: float,
          t_1: float, t_2: float, t: float):
    """Угол тангажа [рад]"""

    if t <= T_V:
        # print(t, np.pi / 2.)
        return np.pi / 2.
    elif T_V <= t and t < t_1:
        # print(t, np.pi / 2. + theta_1 * (t - T_V))
        return np.pi / 2. + theta_1 * (t - T_V)
    elif t_1 <= t and t < t_2:
        return 0
    elif t_2 <= t:
        # print(t, theta_2)
        return theta_2


def FUNCTION_RIGHT_SIDE(u: list, t: float, theta_1: float, theta_2: float,
                        t_1: float, t_2: float, P: Callable, W_EFF: float):
    """Вектор правой части системы"""
    dVx_dt = P(t, t_1, t_2) / u[4] * np.cos(THETA(theta_1,
                                                  theta_2, t_1, t_2, t)) + G_X(u[2], u[3])
    dVy_dt = P(t, t_1, t_2) / u[4] * np.sin(THETA(theta_1,
                                                  theta_2, t_1, t_2, t)) + G_Y(u[2], u[3])
    dx_dt = u[0]
    dy_dt = u[1]
    dm_dt = - P(t, t_1, t_2) / W_EFF

    return dVx_dt, dVy_dt, dx_dt, dy_dt, dm_dt


def RADIUS_VECTOR_CONDITION(x: float, y: float):
    """Проверка условия на радиус вектор"""
    return np.sqrt(x ** 2 + (y + R_M) ** 2) - R_M - H


def SPEED_CONDITION(V_x: float, V_y: float):
    """Проверка условия на скорость"""
    return V_x ** 2 + V_y ** 2 - MU_M / (H + R_M)


def ANGLE_CONDITION(V_x: float, V_y: float, x: float, y: float):
    """
    Проверка условия на угол
    наклона траектории к местному горизонту
    """
    r = pow(x ** 2 + (y + R_M) ** 2, 0.5)
    V = pow(V_x ** 2 + V_y ** 2, 0.5)
    if V < 1e-7:
        return 0
    rV = x * V_x + (y + R_M) * V_y
    return np.arcsin(rV / (r * V))


def THETA_C(V_x: float, V_y: float):
    return np.atan2(V_y, V_x)


def PHI(x: float, y: float):
    if ((R_M ** 2 + R_M * y) /
            (R_M * pow((x ** 2 + (R_M + y) ** 2), 0.5))) >= 1:
        return np.arccos(1)
    return np.arccos((R_M ** 2 + R_M * y) / (R_M * pow((x ** 2 + (R_M + y) ** 2), 0.5)))


def ALPHA(V_x: float, V_y: float, theta: float):
    V = pow((V_x ** 2 + V_y ** 2), 0.5)
    if V < 1e-7:
        return np.arccos(1)
    return np.arccos(V_x / V * np.cos(theta) + V_y / V * np.sin(theta)) / pow((V_x / V) ** 2 + (V_y / V) ** 2, 0.5)
