"""
Compute I = ∫₂³ (x³+1)/(x³-1) dx via scipy.integrate.quad.
"""

import numpy as np
from scipy import integrate


def integrand(x):
    return (x**3 + 1) / (x**3 - 1)


result, error = integrate.quad(integrand, 2, 3)
print("I = ∫₂³ (x³+1)/(x³-1) dx")
print(f"  scipy.quad : {result:.15f}  (est. error {error:.2e})")

result2, error2 = integrate.quad(lambda x: 1 / np.tanh(np.log(x ** (3 / 2))), 2, 3)
print("\nJ = ∫₂³ coth(ln(x^(3/2))) dx")
print(f"  scipy.quad : {result2:.15f}  (est. error {error2:.2e})")
print(f"\n  Discrepancy: {abs(result - result2):.2e}")

result3, error3 = integrate.quad(lambda x: 1 + 2 / (x**3 - 1), 2, 3)
print("\nK = ∫₂³ [1 + 2/(x³-1)] dx")
print(f"  scipy.quad : {result3:.15f}  (est. error {error3:.2e})")
print(f"\n  Discrepancy: {abs(result - result3):.2e}")

# Partial fraction decomposition: 1/(x³-1) = (1/3)/(x-1) - (1/3)(x+2)/(x²+x+1)
# So: 1 + 2/(x³-1) = 1 + (2/3)/(x-1) - (2/3)(x+2)/(x²+x+1)
pf = lambda x: 1 + (2 / 3) / (x - 1) - (2 / 3) * (x + 2) / (x**2 + x + 1)
result4, error4 = integrate.quad(pf, 2, 3)
print("\nL = ∫₂³ [1 + (2/3)/(x-1) - (2/3)(x+2)/(x²+x+1)] dx")
print(f"  scipy.quad : {result4:.15f}  (est. error {error4:.2e})")
print(f"\n  Discrepancy: {abs(result - result4):.2e}")

# --- Exact closed form via SymPy ---
from sympy import symbols, integrate as sym_integrate, simplify, Rational

x = symbols("x")
f = 1 + Rational(2, 3) / (x - 1) - Rational(2, 3) * (x + 2) / (x**2 + x + 1)
antideriv = sym_integrate(f, x)
print("\n--- Exact antiderivative (SymPy) ---")
print(f"  F(x) = {antideriv}")

exact = antideriv.subs(x, 3) - antideriv.subs(x, 2)
exact_simplified = simplify(exact)
print("\n--- Exact value ---")
print(f"  I = {exact_simplified}")
print(f"  I = {float(exact_simplified):.15f}")
print(f"\n  Discrepancy vs quad: {abs(float(exact_simplified) - result):.2e}")
