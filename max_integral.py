"""
Numerical integration of max(x, x^2, x^3) over [-2, 2].

Layer 1: Direct scipy.integrate.quad numerical reference.
"""

from scipy.integrate import quad
from sympy import integrate, symbols


def integrand(x):
    return max(x, x**2, x**3)


result, error = quad(integrand, -2, 2)
print(f"Numerical result: {result:.15f}  err={error:.2e}")
print()

# Verify piecewise decomposition numerically
r1, _ = quad(lambda x: x**2, -2, 0)
r2, _ = quad(lambda x: x, 0, 1)
r3, _ = quad(lambda x: x**3, 1, 2)
piecewise = r1 + r2 + r3
print(f"Piecewise sum:    {piecewise:.15f}")
print(f"  x^2 on [-2,0]: {r1:.15f}")
print(f"  x   on [0,1]:  {r2:.15f}")
print(f"  x^3 on [1,2]:  {r3:.15f}")
print()

# Exact closed form via SymPy
x = symbols("x")
exact = (
    integrate(x**2, (x, -2, 0)) + integrate(x, (x, 0, 1)) + integrate(x**3, (x, 1, 2))
)
print(f"Exact (SymPy):    {exact} = {float(exact):.15f}")
print(f"Difference:       {abs(result - float(exact)):.2e}")
