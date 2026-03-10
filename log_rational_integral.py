"""
Numerical integration of f(x) = ((ln(x) - 1) / (1 + ln²(x)))² over [1, e].

Substitution u = ln(x) gives the equivalent form:
  int_0^1 ((u-1)/(1+u^2))^2 * e^u du

Partial fractions:
  ((u-1)/(1+u^2))^2 = 1/(u^2+1) - 2u/(u^2+1)^2

IBP on the second piece using d/du(1/(u^2+1)) = -2u/(u^2+1)^2:
  2*int_0^1 u*e^u/(u^2+1)^2 du = -[e^u/(u^2+1)]_0^1 + int_0^1 e^u/(u^2+1) du

The int_0^1 e^u/(u^2+1) terms cancel, leaving:
  I = [e^u/(u^2+1)]_0^1 = e/2 - 1
"""

import numpy as np
from scipy.integrate import quad
from sympy import symbols, apart, latex

# --- Layer 1: numerical reference ---


def integrand(x):
    t = np.log(x)
    return ((t - 1) / (1 + t**2)) ** 2


def integrand_sub(u):
    return ((u - 1) / (1 + u**2)) ** 2 * np.exp(u)


result, error = quad(integrand, 1, np.e)
result2, error2 = quad(integrand_sub, 0, 1)
print(f"Layer 1 (original):    {result:.15f}  err={error:.2e}")
print(f"Layer 1 (substituted): {result2:.15f}  err={error2:.2e}")
print(f"Difference:            {abs(result - result2):.2e}")
print()

# --- Partial fraction expansion ---

u = symbols("u")
pf = apart(((u - 1) / (1 + u**2)) ** 2, u)
print("Partial fractions of ((u-1)/(1+u^2))^2:")
print(f"  {pf}")
print(f"  LaTeX: {latex(pf)}")
print()

# --- Closed form via IBP ---

boundary = np.e / 2 - 1
print(f"Closed form e/2 - 1 = {boundary:.15f}")
print(f"Numerical result     = {result:.15f}")
print(f"Difference           = {abs(boundary - result):.2e}")
