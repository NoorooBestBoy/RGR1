import matplotlib.pyplot as plt
import numpy as np
import sympy as sp

rho = 8*10**3
F = 8*10**-3
E = 2.2*10**11
J = 7*10**-6
L = 3
C0 = 1


def calculate_coefficients(k_val):
    C1, C2, C3, C4, k, x = sp.symbols("C1, C2, C3, C4, k, x")
    A = C1 * sp.sin(k * x) + C2 * sp.cos(k * x) + C3 * sp.sinh(k * x) + C4 * sp.cosh(k * x)

    A1 = A.evalf(subs={k: k_val, C4: C0, x: 0})
    A2 = A.diff(x, x).evalf(subs={k: k_val, C4: C0, x: 0})
    A3 = A.evalf(subs={k: k_val, C4: C0, x: L})
    A4 = A.diff(x, x).evalf(subs={k: k_val, C4: C0, x: L})

    res = list(sp.linsolve([A1, A3, A4], (C1, C2, C3)))
    return res[0][0], res[0][1], res[0][2], C0


def calculate_p_roots():
    C1, C2, C3, C4, k, x = sp.symbols("C1, C2, C3, C4, k, x")
    A = C1 * sp.sin(k * x) + C2 * sp.cos(k * x) + C3 * sp.sinh(k * x) + C4 * sp.cosh(k * x)

    A1 = A.evalf(subs={x: 0})
    A2 = A.diff(x, x).evalf(subs={x: 0})
    A3 = A.evalf(subs={x: L})
    A4 = A.diff(x, x).evalf(subs={x: L})
    M = sp.linear_eq_to_matrix([A1, A2, A3, A4], (C1, C2, C3, C4))[0]
    det_M = M.det()

    h = 0.1
    steps_number = 15500
    p0 = 0.1

    p = p0
    x, y = [], []
    roots = []
    for i in range(steps_number):
        k_val = fk(p)
        d = det_M.evalf(subs={k: k_val})

        x.append(p)
        y.append(d)

        if len(y) > 1:
            if y[-2] < 0 < y[-1] or y[-2] > 0 > y[-1]:
                roots.append(x[-1])

        p += h

    print(roots)

    roots_number = 3

    return roots[:roots_number]


def u(x, p, C1, C2, C3, C4):
    k = fk(p)
    res = C1*np.sin(k*x)+C2*np.cos(k*x)+C3*np.sinh(k*x)+C4*np.cosh(k*x)
    return res


def fk(p):
    res = np.power((rho*F*p**2)/(E*J), 1/4)
    return res


def calculate_u():
    roots = calculate_p_roots()
    h = 0.001
    steps_number = int(L / h)

    xs = [i*h for i in range(steps_number)]

    fig, axes = plt.subplots(3, 1) # выводит 3 граф на одном

    for i in range(len(roots)):
        k = fk(roots[i])
        C1, C2, C3, C4 = calculate_coefficients(k)
        ys = [u(x, roots[i], C1, C2, C3, C4) for x in xs]
        axes[i].plot(xs, ys)
        axes[i].grid()
    plt.show()

calculate_u()