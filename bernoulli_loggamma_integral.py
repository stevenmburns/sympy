"""
Numerical integration of B_2(x) * ln(Gamma(x)) over [0, 1],
where B_2(x) = x^2 - x + 1/6 is the second Bernoulli polynomial.

Layer 1: Direct scipy.integrate.quad numerical reference.
"""

import math

import numpy as np
from scipy.integrate import quad
from scipy.special import gammaln, zeta


def integrand(x):
    B2 = x**2 - x + 1 / 6
    return B2 * gammaln(x)


result, error = quad(integrand, 0, 1, limit=100)
print(f"Numerical result:  {result:.15f}  err={error:.2e}")

closed_form = zeta(3) / (4 * math.pi**2)
print(f"zeta(3)/(4*pi^2): {closed_form:.15f}")
print(f"Difference:        {abs(result - closed_form):.2e}")
print()


# --- Step 1: reflection formula ---
# B_2(1-x) = B_2(x), so reflection gives 2I = -int_0^1 B_2(x) ln(sin(pi*x)) dx
def integrand_step1(x):
    B2 = x**2 - x + 1 / 6
    return -0.5 * B2 * np.log(np.sin(np.pi * x))


result_step1, error_step1 = quad(integrand_step1, 0, 1, limit=100)
print(
    f"Step 1 (-1/2 * int B2 ln(sin(pi x))): {result_step1:.15f}  err={error_step1:.2e}"
)
print(f"Matches original:                      {abs(result - result_step1):.2e}")
