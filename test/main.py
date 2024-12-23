import math
import numpy as np

# Constants
MU_M = 4903.0 * 10**9
R = 1738000.0
TV = 14
P1 = 9980.0  # or 8300.0 if needed

# Define the Part class (equivalent to C++ struct)


class Part:
    def __init__(self, Vx=0.0, Vy=0.0, x=0.0, y=0.0, m=0.0, t=0.0):
        self.Vx = Vx
        self.Vy = Vy
        self.x = x
        self.y = y
        self.m = m
        self.t = t

# Helper function to calculate Teta


def Teta_func(x, y, Vx, Vy):
    rV = x * Vx + (y + R) * Vy
    r = math.sqrt(x**2 + (y + R)**2)
    V = math.sqrt(Vx**2 + Vy**2)
    return math.asin(rV / (r * V))

# Helper function to calculate RF


def RF_func(x, y, h_isl):
    r = math.sqrt(x**2 + (y + R)**2)
    return r

# Helper function to calculate Teta_C


def Teta_C(Vx, Vy):
    return math.atan2(Vy, Vx)

# Function RP (equivalent to C++ RP function)


def RP(Left, delta_dTetta_dt, delta_Tetta2, flag, dTeta_dt, Teta2, Wist):
    Res = Part()

    delta = 1 if flag in [0, 1, 2] else 0
    Teta = math.pi / 2  # Default Teta
    if flag == 1:
        Teta = math.pi / 2 + dTeta_dt * (Left.t - TV)
    elif flag == 2:
        Teta = Teta2

    beta = P1 * delta / Wist
    gx = -MU_M * Left.x / (math.sqrt(Left.x**2 + (Left.y + R)**2))**3
    gy = -MU_M * (Left.y + R) / (math.sqrt(Left.x**2 + (Left.y + R)**2))**3

    Res.Vx = delta * P1 / Left.m * math.cos(Teta) + gx
    Res.Vy = delta * P1 / Left.m * math.sin(Teta) + gy
    Res.x = Left.Vx * 1
    Res.y = Left.Vy * 1
    Res.m = -beta
    Res.t = 1.0

    return Res

# Function to sum two Parts


def PLUS(P1, P2):
    P3 = Part()
    P3.Vx = P1.Vx + P2.Vx
    P3.Vy = P1.Vy + P2.Vy
    P3.x = P1.x + P2.x
    P3.y = P1.y + P2.y
    P3.m = P1.m + P2.m
    P3.t = P1.t + P2.t
    return P3

# Function to multiply a Part by a constant


def UMN(P1, c):
    P = Part()
    P.Vx = P1.Vx * c
    P.Vy = P1.Vy * c
    P.x = P1.x * c
    P.y = P1.y * c
    P.m = P1.m * c
    P.t = P1.t * c
    return P

# RK4_step function in Python


def RK4_step(Left, dt, delta_Tetta_dt, delta_Tetta2, flag, dTeta_dt, Teta2, Wist):
    K = [None] * 4
    Left0 = Left
    delta = 0
    Teta = 0

    for i in range(4):
        K[i] = RP(Left0, delta_Tetta_dt, delta_Tetta2,
                  flag, dTeta_dt, Teta2, Wist)
        if i <= 1:
            delta = dt / 2.0
        if i == 2:
            delta = dt
        if i != 3:
            Left0 = PLUS(Left, UMN(K[i], delta))

    Left = PLUS(Left, UMN(K[0], dt * 1.0 / 6.0))
    Left = PLUS(Left, UMN(K[1], dt * 2.0 / 6.0))
    Left = PLUS(Left, UMN(K[2], dt * 2.0 / 6.0))
    Left = PLUS(Left, UMN(K[3], dt * 1.0 / 6.0))

    return Left, Teta

# RK4 function in Python


def RK4(Left, dt, delta_Tetta_dt, delta_Tetta2, flag, dTeta_dt, Teta2, delta_r, delta_Teta, delta_m, t1, t2, h_isl1_1, Wist, m_k, m_pg):
    tv = 14
    flag_1 = True
    flag_2 = True
    dt1, dt2 = 0.0, 0.0
    R = 1738000.0
    mu_m = 4903.0 * 10**9
    R_isl1_1 = h_isl1_1 + R
    Vk = math.sqrt(mu_m / R_isl1_1)
    Tang = math.pi / 2  # Тангент угла
    Buf = Part()

    r = RF_func(Left.x, Left.y, h_isl1_1)
    modV = math.sqrt(Left.Vx**2 + Left.Vy**2)
    TetaC = Teta_C(Left.Vx, Left.Vy)
    Teta = Teta_func(Left.x, Left.y, Left.Vx, Left.Vy)
    h = r - R
    alpha = Tang - TetaC
    Fi = math.acos((Left.y + R) / r)

    while math.sqrt(Left.Vx**2 + Left.Vy**2) < Vk:
        st = str(int(Left.t))

        if Left.t < 99:
            if len(st) < 3 or st == "14":
                pass  # Можно реализовать вывод данных

        if Left.t >= 99 and len(st) < 4:
            pass

        if Left.t >= 999 and len(st) < 5:
            pass

        Buf = Left
        RK4_step(Left, dt, delta_Tetta_dt, delta_Tetta2,
                 flag, dTeta_dt, Teta2, Wist, Tang)

        r = RF_func(Left.x, Left.y, h_isl1_1)
        modV = math.sqrt(Left.Vx**2 + Left.Vy**2)
        TetaC = Teta_C(Left.Vx, Left.Vy)
        Teta = Teta_func(Left.x, Left.y, Left.Vx, Left.Vy)
        h = r - R
        alpha = Tang - TetaC
        Fi = math.acos((Left.y + R) / r)

        if Left.t > tv and Left.t < t1:
            flag = 1

        if flag_1 and Left.t > t1 and Left.t < t2:
            Left = Buf
            dt1 = t1 - Left.t
            RK4_step(Left, dt1, delta_Tetta_dt, delta_Tetta2,
                     flag, dTeta_dt, Teta2, Wist, Tang)

            r = RF_func(Left.x, Left.y, h_isl1_1)
            modV = math.sqrt(Left.Vx**2 + Left.Vy**2)
            TetaC = Teta_C(Left.Vx, Left.Vy)
            Teta = Teta_func(Left.x, Left.y, Left.Vx, Left.Vy)
            h = r - R
            alpha = Tang - TetaC
            Fi = math.acos((Left.y + R) / r)

            dt1 = dt - dt1
            flag = 3
            RK4_step(Left, dt1, delta_Tetta_dt, delta_Tetta2,
                     flag, dTeta_dt, Teta2, Wist, Tang)

            r = RF_func(Left.x, Left.y, h_isl1_1)
            modV = math.sqrt(Left.Vx**2 + Left.Vy**2)
            TetaC = Teta_C(Left.Vx, Left.Vy)
            Teta = Teta_func(Left.x, Left.y, Left.Vx, Left.Vy)
            h = r - R
            alpha = Tang - TetaC
            Fi = math.acos((Left.y + R) / r)

            flag_1 = False

        if flag_2 and Left.t > t2:
            Left = Buf
            dt2 = t2 - Left.t
            RK4_step(Left, dt2, delta_Tetta_dt, delta_Tetta2,
                     flag, dTeta_dt, Teta2, Wist, Tang)

            r = RF_func(Left.x, Left.y, h_isl1_1)
            modV = math.sqrt(Left.Vx**2 + Left.Vy**2)
            TetaC = Teta_C(Left.Vx, Left.Vy)
            Teta = Teta_func(Left.x, Left.y, Left.Vx, Left.Vy)
            h = r - R
            alpha = Tang - TetaC
            Fi = math.acos((Left.y + R) / r)

            dt2 = dt - dt2
            flag = 2
            RK4_step(Left, dt2, delta_Tetta_dt, delta_Tetta2,
                     flag, dTeta_dt, Teta2, Wist, Tang)

            r = RF_func(Left.x, Left.y, h_isl1_1)
            modV = math.sqrt(Left.Vx**2 + Left.Vy**2)
            TetaC = Teta_C(Left.Vx, Left.Vy)
            Teta = Teta_func(Left.x, Left.y, Left.Vx, Left.Vy)
            h = r - R
            alpha = Tang - TetaC
            Fi = math.acos((Left.y + R) / r)

            flag_2 = False

    if abs(modV - Vk) > 1.0E-06:
        while abs(modV - Vk) > 1.0E-06:
            dt = -dt
            RK4_step(Left, dt, delta_Tetta_dt, delta_Tetta2,
                     flag, dTeta_dt, Teta2, Wist, Tang)

            r = RF_func(Left.x, Left.y, h_isl1_1)
            modV = math.sqrt(Left.Vx**2 + Left.Vy**2)
            TetaC = Teta_C(Left.Vx, Left.Vy)
            Teta = Teta_func(Left.x, Left.y, Left.Vx, Left.Vy)
            h = r - R
            alpha = Tang - TetaC
            Fi = math.acos((Left.y + R) / r)

            dt = -dt / 5.0
            modV = math.sqrt(Left.Vx**2 + Left.Vy**2)

            while modV < Vk:
                RK4_step(Left, dt, delta_Tetta_dt, delta_Tetta2,
                         flag, dTeta_dt, Teta2, Wist, Tang)

                r = RF_func(Left.x, Left.y, h_isl1_1)
                TetaC = Teta_C(Left.Vx, Left.Vy)
                Teta = Teta_func(Left.x, Left.y, Left.Vx, Left.Vy)
                h = r - R
                alpha = Tang - TetaC
                Fi = math.acos((Left.y + R) / r)

                modV = math.sqrt(Left.Vx**2 + Left.Vy**2)

    delta_r[0] = RF_func(Left.x, Left.y, h_isl1_1)
    delta_Teta[0] = Teta_func(Left.x, Left.y, Left.Vx, Left.Vy)
    delta_m[0] = 3000.0 - Left.m
    m_pg[0] = Left.m - m_k


Left = Part(Vx=5000, Vy=3000, x=1000, y=2000, m=500, t=0.0)
RK4(Left, dt=0.1, delta_Tetta_dt=0.05, delta_Tetta2=0.02, flag=1, dTeta_dt=0.03, Teta2=0.5, delta_r=None,
    delta_Teta=None, delta_m=None, t1=50, t2=100, h_isl1_1=300, Wist=0.01, m_k=0.01, m_pg=0.02)
