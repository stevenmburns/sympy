import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

x = np.linspace(0, 2, 1000)

phi = (1 + np.sqrt(5)) / 2

endpoints = [(0, 1 / phi), (1 / phi, 1), (1, phi), (phi, 2)]

x0 = np.linspace(*endpoints[0], 500)
x1 = np.linspace(*endpoints[1], 500)
x2 = np.linspace(*endpoints[2], 500)
x3 = np.linspace(*endpoints[3], 500)


def f(x, r, s):
    return r * (s * (x**3 - x) - x**2)


def F(x, r, s):
    return r * (s * (x**4 / 4 - x**2 / 2) - x**3 / 3)


def g(x):
    return np.abs(np.abs(x**3 - x) - x**2)


result, error = integrate.quad(g, 0, 2)
print(f"∫₀² g(x) dx ≈ {result:.10f}  (estimated error: {error:.2e})")

y = g(x)

params = [(1, -1), (-1, -1), (-1, 1), (1, 1)]

integrals = [F(b, *p) - F(a, *p) for (a, b), p in zip(endpoints, params)]
for i, ((a, b), v) in enumerate(zip(endpoints, integrals)):
    print(f"∫_{a:.6f}^{b:.6f} f(x, r={params[i][0]}, s={params[i][1]}) dx = {v:.10f}")

phi_sym = (1 + sqrt(5)) / 2
sym_endpoints = [
    (Rational(0), 1 / phi_sym),
    (1 / phi_sym, Rational(1)),
    (Rational(1), phi_sym),
    (phi_sym, Rational(2)),
]


def F_sym(x, r, s):
    return r * (s * (x**4 / 4 - x**2 / 2) - x**3 / 3)


print("\nExact symbolic integrals:")
for i, ((a, b), p) in enumerate(zip(sym_endpoints, params)):
    exact = simplify(F_sym(b, *p) - F_sym(a, *p))
    print(f"  segment {i}: {exact} ≈ {float(exact):.10f}")

y0 = f(x0, *params[0])
y1 = f(x1, *params[1])
y2 = f(x2, *params[2])
y3 = f(x3, *params[3])

plt.plot(x, y, label=r"$|\ |x^3-x|\ -x^2\ |$")
plt.plot(x0, y0, label=r"$f(x,\ r{=}1,\ s{=}{-1})$")
plt.plot(x1, y1, label=r"$f(x,\ r{=}{-1},\ s{=}{-1})$")
plt.plot(x2, y2, label=r"$f(x,\ r{=}{-1},\ s{=}1)$")
plt.plot(x3, y3, label=r"$f(x,\ r{=}1,\ s{=}1)$")
plt.legend()
plt.xlabel("x")
plt.ylabel("y")
plt.title(r"$|\ |x^3 - x|\ - x^2\ |$")
plt.grid(True)
plt.tight_layout()
# plt.savefig("plot.png", dpi=150)
plt.show()
