"""
Compute ∫₀¹ arctan(x)/x · ln(1+x) dx via a digamma series.

Derivation:
  ln(1+x) = Σ_{k=1}^∞ (-1)^{k-1} x^k / k

  So  I = Σ_{k=1}^∞ (-1)^{k-1}/k · ∫₀¹ x^{k-1} arctan(x) dx

Each integral I_k = ∫₀¹ x^{k-1} arctan(x) dx is evaluated by:
  1. IBP:  I_k = π/(4k) - (1/k) ∫₀¹ x^k/(1+x²) dx
  2. Sub t=x²: ∫₀¹ x^k/(1+x²) dx = (1/2) ∫₀¹ t^{(k-1)/2}/(1+t) dt
  3. Formula ∫₀¹ t^{s-1}/(1+t) dt = (1/2)[ψ((s+1)/2) - ψ(s/2)], s=(k+1)/2:
           = (1/4)[ψ((k+3)/4) - ψ((k+1)/4)]

  => I_k = (1/(4k)) [π + ψ((k+1)/4) - ψ((k+3)/4)]

  => I = Σ_{k=1}^∞ (-1)^{k-1} / (4k²) · [π + ψ((k+1)/4) - ψ((k+3)/4)]
"""

import numpy as np
from scipy import integrate as sci_integrate
from scipy.special import digamma


# --- Reference via scipy.quad ---
result_ref, _ = sci_integrate.quad(
    lambda t: np.arctan(t) / t * np.log(1 + t) if t > 0 else 0.0, 0, 1
)
print("∫₀¹ arctan(x)/x · ln(1+x) dx")
print(f"  scipy.quad reference : {result_ref:.15f}\n")


def digamma_term(k):
    """I_k = (1/(4k)) [π + ψ((k+1)/4) - ψ((k+3)/4)]"""
    return (np.pi + digamma((k + 1) / 4) - digamma((k + 3) / 4)) / (4 * k)


def series_partial(N):
    """Σ_{k=1}^N (-1)^{k-1}/k · I_k"""
    return sum(((-1) ** (k - 1) / k) * digamma_term(k) for k in range(1, N + 1))


# Verify I_1 = π/4 - ln2/2 ≈ 0.43882 analytically
I1_formula = digamma_term(1)
I1_exact = np.pi / 4 - np.log(2) / 2
print("Verification I_1 = ∫₀¹ arctan(x) dx:")
print(f"  Formula  : {I1_formula:.15f}")
print(f"  Exact π/4 - ln2/2 : {I1_exact:.15f}\n")

print("Digamma series convergence:")
print(f"  {'N':>6}  {'partial sum':>20}  {'error':>12}")
for N in [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000]:
    s = series_partial(N)
    print(f"  {N:>6}  {s:>20.15f}  {abs(s - result_ref):>12.2e}")
