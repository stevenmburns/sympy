"""
Numerical integration of the infinite nested radical over [0, 1]:

  K(x) = sqrt(x^2 + x + sqrt(x^2 + x + sqrt(x^2 + x + ...)))

K satisfies K^2 = x^2 + x + K, whose positive root is K = x + 1.
So the integral equals int_0^1 (x+1) dx = 3/2.

Layer 1: Direct numerical reference — iterate the nested radical to convergence.
"""

import numpy as np
from scipy.integrate import quad


def nested_radical(x, depth=100):
    """Evaluate K(x) by iterating K <- sqrt(x^2 + x + K) from K=0."""
    K = 0.0
    for _ in range(depth):
        K = np.sqrt(x**2 + x + K)
    return K


result, error = quad(nested_radical, 0, 1)
print(f"Numerical result (nested radical): {result:.15f}  err={error:.2e}")
print(f"Exact 3/2:                         {3 / 2:.15f}")
print(f"Difference:                        {abs(result - 1.5):.2e}")
