import numpy as np
from runge_kutta import runge_kutta_4


def test_runge_kutta_exponential_growth():
    """
    Тест 1: Проверка на экспоненциальное уравнение 
    dy/dt = y с аналитическим решением y(t) = e^t
    """
    print("[INFO] Make test test_exponential_growth.")

    def exponential_rhs(u, t):
        return [u[0]]  # правая часть уравнения, возвращает производную y

    # Начальные условия
    u_0 = [1.0]  # y(0) = 1
    t_0 = 0.0  # Начальное время
    t_k = 1.0  # Конечное время
    tau = 1e-2  # Шаг интегрирования

    # Численное решение с помощью метода Рунге-Кутты
    solution = runge_kutta_4(u_0, t_0, t_k, tau, exponential_rhs)

    # Получаем последнее численное значение y(t) и сравниваем с аналитическим решением
    y_numeric = solution[-1][0]  # Численное значение в t_k
    y_analytic = np.exp(t_k)  # Аналитическое значение в t_k

    # Допустимая погрешность
    tol = 1e-2
    assert abs(y_numeric - y_analytic) < tol, f"[ERROR] Test failed: y_numeric = {
        y_numeric}, y_analytic = {y_analytic}"


def test_runge_kutta_harmonic_oscillator():
    """
    Тест 2: Проверка на простое 
    гармоническое уравнение
    """
    print("[INFO] Make test test_harmonic_oscillator.")
    # Задаем правую часть уравнения: d²x/dt² = -x или в виде системы уравнений

    def harmonic_rhs(u, t):
        x, v = u  # u = [x, v], где x - положение, v - скорость
        return [v, -x]  # Возвращаем [dx/dt, dv/dt]

    # Начальные условия
    u_0 = [1.0, 0.0]  # x(0) = 1, v(0) = 0
    t_0 = 0.0  # Начальное время
    # Время, равное периоду (2*pi для гармонического осциллятора с w=1)
    t_k = 2 * np.pi
    tau = 1e-2  # Шаг интегрирования

    # Численное решение с помощью метода Рунге-Кутты
    solution = runge_kutta_4(u_0, t_0, t_k, tau, harmonic_rhs)

    # Получаем последнее численное значение x(t) и сравниваем с аналитическим решением
    x_numeric = solution[-1][0]  # Численное значение x в t_k
    x_analytic = 1.0  # Аналитическое значение x(t) = cos(t) при t = 2*pi

    # Допустимая погрешность
    tol = 1e-2
    assert abs(x_numeric - x_analytic) < tol, f"[ERROR] Test failed: x_numeric = {
        x_numeric}, x_analytic = {x_analytic}"


if __name__ == "__main__":
    # Запуск тестов
    test_runge_kutta_exponential_growth()
    test_runge_kutta_harmonic_oscillator()
    print("[INFO] All tests passed.")
