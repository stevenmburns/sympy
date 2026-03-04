"""
Compute ∫₀¹ arctan(x)/x · ln(1+x) dx via two digamma series.

Series 1 — expand ln(1+x):
  I_k = ∫₀¹ x^{k-1} arctan(x) dx
      = (1/(4k)) [π + ψ((k+1)/4) - ψ((k+3)/4)]
  I = Σ_{k=1}^∞ (-1)^{k-1} / (4k²) · [π + ψ((k+1)/4) - ψ((k+3)/4)]

Series 2 — expand arctan(x)/x:
  J_n = ∫₀¹ x^{2n} ln(1+x) dx
      = (1/(2n+1)) [ln2 - (1/2)(ψ(n+3/2) - ψ(n+1))]
  I = Σ_{n=0}^∞ (-1)^n / (2n+1)² · [ln2 - (1/2)(ψ(n+3/2) - ψ(n+1))]

Both series are accelerated with iterated Aitken Δ² applied to partial sums.
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


# --- Series terms ---


def series1_term(k):
    """(-1)^{k-1} / (4k²) · [π + ψ((k+1)/4) - ψ((k+3)/4)]"""
    bracket = np.pi + digamma((k + 1) / 4) - digamma((k + 3) / 4)
    return ((-1) ** (k - 1)) / (4 * k**2) * bracket


def series2_term(n):
    """(-1)^n / (2n+1)² · [ln2 - (1/2)(ψ(n+3/2) - ψ(n+1))]"""
    bracket = np.log(2) - 0.5 * (digamma(n + 1.5) - digamma(n + 1))
    return ((-1) ** n) / (2 * n + 1) ** 2 * bracket


# --- Aitken Δ² acceleration ---
# Given partial sums S_n, S_{n+1}, S_{n+2}:
#   Aitken(S_n) = S_n - (S_{n+1} - S_n)² / (S_{n+2} - 2·S_{n+1} + S_n)


def aitken_accelerate(partial_sums):
    """One pass of Aitken Δ² over a sequence of partial sums."""
    acc = []
    for i in range(len(partial_sums) - 2):
        s0, s1, s2 = partial_sums[i], partial_sums[i + 1], partial_sums[i + 2]
        denom = s2 - 2 * s1 + s0
        if abs(denom) < 1e-30:
            acc.append(s0)
        else:
            acc.append(s0 - (s1 - s0) ** 2 / denom)
    return acc


def run_series(label, terms):
    """Print convergence table and iterated Aitken acceleration for a term list."""
    # Partial sums
    psums = []
    s = 0.0
    for t in terms:
        s += t
        psums.append(s)

    N_vals = [1, 2, 5, 10, 20, 50, len(terms)]
    print(f"{label} — raw convergence:")
    print(f"  {'N':>6}  {'partial sum':>20}  {'error':>12}")
    for N in N_vals:
        print(
            f"  {N:>6}  {psums[N - 1]:>20.15f}  {abs(psums[N - 1] - result_ref):>12.2e}"
        )

    print(f"\n{label} + iterated Aitken Δ²:")
    print(f"  {'pass':>6}  {'last estimate':>20}  {'error':>12}")
    seq = psums[:]
    for p in range(7):
        print(f"  {p:>6}  {seq[-1]:>20.15f}  {abs(seq[-1] - result_ref):>12.2e}")
        seq = aitken_accelerate(seq)
        if len(seq) < 3:
            break
    print()


# --- Verify first terms analytically ---
# series1_term(1) = (1/4)(π - 2ln2) = π/4 - ln2/2 = ∫₀¹ arctan(x) dx
print("Verification: series1_term(1) = ∫₀¹ arctan(x) dx = π/4 - ln2/2:")
print(f"  series1_term(1) : {series1_term(1):.15f}")
print(f"  exact π/4-ln2/2 : {np.pi / 4 - np.log(2) / 2:.15f}\n")

# series2_term(0) = ln2 - (1/2)(ψ(3/2) - ψ(1)) = 2ln2 - 1 = ∫₀¹ ln(1+x) dx
print("Verification: series2_term(0) = ∫₀¹ ln(1+x) dx = 2ln2 - 1:")
print(f"  series2_term(0) : {series2_term(0):.15f}")
print(f"  exact 2·ln2 - 1 : {2 * np.log(2) - 1:.15f}\n")

MAX_N = 50
run_series("Series 1 (expand ln(1+x))", [series1_term(k) for k in range(1, MAX_N + 1)])
run_series("Series 2 (expand arctan(x))", [series2_term(n) for n in range(MAX_N)])
