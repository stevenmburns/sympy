"""
Compute I = ∫₀¹ sin(πx) / (sin(πx)·Γ²(x) + π) dx via scipy.integrate.quad.
"""

import numpy as np
from scipy import integrate
from scipy.special import gamma


def integrand(x):
    s = np.sin(np.pi * x)
    g = gamma(x)
    return s / (s * g**2 + np.pi)


result, error = integrate.quad(integrand, 0, 1)
print("I = ∫₀¹ sin(πx) / (sin(πx)·Γ²(x) + π) dx")
print(f"  scipy.quad : {result:.15f}  (est. error {error:.2e})")


# Using Γ(x)Γ(1-x) = π/sin(πx), the integrand becomes 1/(Γ²(x) + Γ(x)Γ(1-x))
def integrand2(x):
    g = gamma(x)
    g1 = gamma(1 - x)
    return 1 / (g**2 + g * g1)


result2, error2 = integrate.quad(integrand2, 0, 1)
print("\nJ = ∫₀¹ 1/(Γ²(x) + Γ(x)Γ(1-x)) dx")
print(f"  scipy.quad : {result2:.15f}  (est. error {error2:.2e})")
print(f"\n  Discrepancy: {abs(result - result2):.2e}")

result3, error3 = integrate.quad(
    lambda x: 1 / (gamma(x) * (gamma(x) + gamma(1 - x))), 0, 1
)
print("\nK = ∫₀¹ 1/(Γ(x)(Γ(x)+Γ(1-x))) dx")
print(f"  scipy.quad : {result3:.15f}  (est. error {error3:.2e})")
print(f"\n  Discrepancy: {abs(result - result3):.2e}")

# King's rule: x -> 1-x gives 1/(Γ(1-x)(Γ(1-x)+Γ(x)))
result4, error4 = integrate.quad(
    lambda x: 1 / (gamma(1 - x) * (gamma(1 - x) + gamma(x))), 0, 1
)
print("\nL = ∫₀¹ 1/(Γ(1-x)(Γ(1-x)+Γ(x))) dx  [king's rule applied to K]")
print(f"  scipy.quad : {result4:.15f}  (est. error {error4:.2e})")
print(f"\n  Discrepancy: {abs(result - result4):.2e}")

# Average of K and L:
# (K+L)/2 = (1/2) ∫₀¹ (Γ(x)+Γ(1-x)) / (Γ(x)Γ(1-x)(Γ(x)+Γ(1-x))) dx
#          = (1/2) ∫₀¹ 1/(Γ(x)Γ(1-x)) dx
# Reflection formula Γ(x)Γ(1-x) = π/sin(πx) gives 1/(Γ(x)Γ(1-x)) = sin(πx)/π
# So I = (1/2π) ∫₀¹ sin(πx) dx
result5, error5 = integrate.quad(lambda x: np.sin(np.pi * x) / (2 * np.pi), 0, 1)
print(
    "\nM = (1/2π) ∫₀¹ sin(πx) dx  [after averaging K and L, applying reflection formula]"
)
print(f"  scipy.quad : {result5:.15f}  (est. error {error5:.2e})")
print(f"\n  Discrepancy: {abs(result - result5):.2e}")

# --- Exact answer via SymPy ---
from sympy import symbols, sin, pi, integrate as sym_integrate, simplify

x = symbols("x")
exact = sym_integrate(sin(pi * x) / (2 * pi), (x, 0, 1))
exact_simplified = simplify(exact)
print("\n--- Exact answer (SymPy) ---")
print(f"  I = {exact_simplified}")
print(f"  I = {float(exact_simplified):.15f}")
print(f"\n  Discrepancy vs quad: {abs(float(exact_simplified) - result):.2e}")
