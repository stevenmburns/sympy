"""
Compute I = ∫₀^∞ e^{-x} sin(x) erf(x) dx via scipy.integrate.quad.

--- Derivation via IBP ---

Let u = erf(x), dv = e^{-x} sin(x) dx.
Then du = (2/√π) e^{-x²} dx, and integrating by parts:

  ∫ e^{-x} sin(x) dx  =  -e^{-x}(sin x + cos x) / 2

Boundary terms: at x=0, erf(0)=0; at x=∞, e^{-x}→0. Both vanish.

  I = ∫₀^∞ (2/√π) e^{-x²} · e^{-x}(sin x + cos x)/2 dx
    = (1/√π) ∫₀^∞ e^{-x²-x}(sin x + cos x) dx

--- Completing the square ---

  -x² - x = -(x + 1/2)² + 1/4,  so  e^{-x²-x} = e^{1/4} e^{-(x+1/2)²}

Substitute t = x + 1/2 (limits shift from 0 → 1/2):

  I = (e^{1/4}/√π) ∫_{1/2}^∞ e^{-t²} [sin(t - 1/2) + cos(t - 1/2)] dt

Since sin θ + cos θ = √2 sin(θ + π/4):

  I = (√2 e^{1/4}/√π) ∫_{1/2}^∞ e^{-t²} sin(t + π/4 - 1/2) dt

--- Closed form via complex erfc ---

Write sin(t + α) = Im[e^{i(t+α)}] with α = π/4 - 1/2:

  ∫_{1/2}^∞ e^{-t²} sin(t + α) dt = Im[ e^{iα} ∫_{1/2}^∞ e^{-t²+it} dt ]

For the inner integral, complete the complex square:
  -t² + it = -(t - i/2)² - 1/4

  ∫_{1/2}^∞ e^{-t²+it} dt = e^{-1/4} ∫_{1/2}^∞ e^{-(t-i/2)²} dt
                           = e^{-1/4} · (√π/2) erfc(1/2 - i/2)

Putting it all together:

  I = (√2 e^{1/4}/√π) · Im[ e^{iα} · e^{-1/4} · (√π/2) erfc(1/2 - i/2) ]
    = (1/√2) Im[ e^{i(π/4 - 1/2)} erfc(1/2 - i/2) ]

where erfc(z) = 1 - erf(z) is the complementary error function at complex argument.

--- Series expansion ---

Starting from the IBP form, split e^{-t²-t} = e^{-t²} · e^{-t} and expand e^{-t}:

  I = (1/√π) Σ_{k=0}^∞ (-1)^k/k! · a_k

where  a_k = ∫₀^∞ t^k e^{-t²}(sin t + cos t) dt.

Since |a_k| ≤ √2 · Γ((k+1)/2)/2 grows only sub-factorially, the term a_k/k!
decays super-exponentially and the series converges extremely fast.
"""

import math

import numpy as np
from scipy import integrate
from scipy.special import erf


def integrand(x):
    return np.exp(-x) * np.sin(x) * erf(x)


# --- Direct numerical integration ---
result, error = integrate.quad(integrand, 0, np.inf)
print("I = ∫₀^∞ e^{-x} sin(x) erf(x) dx")
print(f"  scipy.quad  : {result:.15f}  (est. error {error:.2e})")

# --- IBP reduction ---
# I = (1/√π) ∫₀^∞ e^{-x²-x}(sin x + cos x) dx
result_ibp, error_ibp = integrate.quad(
    lambda x: np.exp(-(x**2) - x) * (np.sin(x) + np.cos(x)), 0, np.inf
)
result_ibp /= np.sqrt(np.pi)
print(
    f"\n  After IBP   : {result_ibp:.15f}  (est. error {error_ibp / np.sqrt(np.pi):.2e})"
)
print("  (1/√π) ∫₀^∞ e^{-x²-x}(sin x + cos x) dx")

# --- Closed form via complex erfc ---
# I = (1/√2) Im[ e^{i(π/4 - 1/2)} erfc(1/2 - i/2) ]
z = 0.5 - 0.5j
erfc_z = 1.0 - erf(z)  # scipy.special.erf supports complex arguments
alpha = np.pi / 4 - 0.5
closed_form = (1.0 / np.sqrt(2)) * np.imag(np.exp(1j * alpha) * erfc_z)

print(f"\n  Closed form : {closed_form:.15f}")
print("  (1/√2) Im[ e^{i(π/4 - 1/2)} erfc(1/2 - i/2) ]")
print(f"  erfc(1/2 - i/2) = {erfc_z.real:.10f} + {erfc_z.imag:.10f}i")

print(f"\n  Discrepancy (quad vs closed form): {abs(result - closed_form):.2e}")

# --- Series expansion ---
# I = (1/√π) Σ_{k=0}^∞ (-1)^k/k! · a_k
# a_k = ∫₀^∞ t^k e^{-t²}(sin t + cos t) dt


def series_term(k):
    a_k, _ = integrate.quad(
        lambda t: t**k * np.exp(-(t**2)) * (np.sin(t) + np.cos(t)), 0, np.inf
    )
    return (-1) ** k / math.factorial(k) / np.sqrt(np.pi) * a_k


def aitken_accelerate(partial_sums):
    acc = []
    for i in range(len(partial_sums) - 2):
        s0, s1, s2 = partial_sums[i], partial_sums[i + 1], partial_sums[i + 2]
        denom = s2 - 2 * s1 + s0
        acc.append(s0 if abs(denom) < 1e-30 else s0 - (s1 - s0) ** 2 / denom)
    return acc


MAX_K = 20
terms = [series_term(k) for k in range(MAX_K + 1)]
partial_sums = list(np.cumsum(terms))

print("\n--- Series: (1/√π) Σ (-1)^k/k! · a_k ---")
print(f"  {'N':>4}  {'partial sum':>20}  {'error':>12}")
for N in [1, 2, 3, 5, 8, 12, MAX_K + 1]:
    ps = partial_sums[N - 1]
    print(f"  {N:>4}  {ps:>20.15f}  {abs(ps - result):>12.2e}")

print("\n  + iterated Aitken Δ²:")
print(f"  {'pass':>5}  {'last estimate':>20}  {'error':>12}")
seq = partial_sums[:]
for p in range(6):
    print(f"  {p:>5}  {seq[-1]:>20.15f}  {abs(seq[-1] - result):>12.2e}")
    seq = aitken_accelerate(seq)
    if len(seq) < 3:
        break
