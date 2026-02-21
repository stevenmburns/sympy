import numpy as np
from scipy import integrate
from sympy import Symbol, summation, Rational, polygamma, pi, sin, simplify, Float


def integrand(x):
    # Numerically stable form of exp(4x/5 - ln2) * sech(x)
    # = exp(4x/5) / (e^x + e^{-x})
    if x >= 0:
        return np.exp(-x / 5) / (1 + np.exp(-2 * x))
    else:
        return np.exp(9 * x / 5) / (1 + np.exp(2 * x))


def asymptotic_pos_inf(x, n_terms=5):
    # x -> +inf: f(x) = exp(-x/5)/(1 + exp(-2x))
    # = sum_{n=0}^{N} (-1)^n * exp(-(2n + 1/5)*x)
    return sum((-1) ** n * np.exp(-(2 * n + 0.2) * x) for n in range(n_terms))


def asymptotic_neg_inf(x, n_terms=5):
    # x -> -inf: f(x) = exp(9x/5)/(1 + exp(2x))
    # = sum_{n=0}^{N} (-1)^n * exp((2n + 9/5)*x)
    return sum((-1) ** n * np.exp((2 * n + 1.8) * x) for n in range(n_terms))


def integrate_series(n_terms):
    # ∫_0^∞ sum_n (-1)^n exp(-(2n+1/5)x) dx = sum_n (-1)^n / (2n + 1/5)
    pos = sum((-1) ** n / (2 * n + 0.2) for n in range(n_terms))
    # ∫_{-∞}^0 sum_n (-1)^n exp((2n+9/5)x) dx = sum_n (-1)^n / (2n + 9/5)
    neg = sum((-1) ** n / (2 * n + 1.8) for n in range(n_terms))
    return pos, neg


result, error = integrate.quad(integrand, -np.inf, np.inf)
print(f"Numerical result: {result}")
print(f"Closed form π/(2 sin(π/10)): {np.pi / (2 * np.sin(np.pi / 10))}")

print("\nAsymptotic expansion vs exact (x -> +inf):")
print(f"{'x':>6}  {'exact':>20}  {'asymptotic':>20}  {'rel error':>12}")
for x in [2, 5, 10, 20]:
    exact = integrand(x)
    approx = asymptotic_pos_inf(x)
    print(f"{x:>6}  {exact:>20.6e}  {approx:>20.6e}  {abs(approx - exact) / exact:>12.2e}")

print("\nAsymptotic expansion vs exact (x -> -inf):")
print(f"{'x':>6}  {'exact':>20}  {'asymptotic':>20}  {'rel error':>12}")
for x in [-2, -5, -10, -20]:
    exact = integrand(x)
    approx = asymptotic_neg_inf(x)
    print(f"{x:>6}  {exact:>20.6e}  {approx:>20.6e}  {abs(approx - exact) / exact:>12.2e}")

print("\nIntegrating asymptotic series term by term:")
print(f"  ∫_0^∞  ~ sum_n (-1)^n / (2n + 1/5)  = 5 sum_n (-1)^n / (10n + 1)")
print(f"  ∫_-∞^0 ~ sum_n (-1)^n / (2n + 9/5)  = 5 sum_n (-1)^n / (10n + 9)")
print(f"\n{'terms':>6}  {'pos half':>16}  {'neg half':>16}  {'total':>16}  {'error':>12}")
for n in [1, 2, 5, 10, 20, 50, 100]:
    pos, neg = integrate_series(n)
    total = pos + neg
    print(f"{n:>6}  {pos:>16.10f}  {neg:>16.10f}  {total:>16.10f}  {abs(total - result):>12.2e}")

# --- Exact closed forms via digamma ---
# sum_{n=0}^inf (-1)^n / (n + a) = (1/2)[psi((a+1)/2) - psi(a/2)]
# so sum_{n=0}^inf (-1)^n / (10n + r) = (1/10) * (1/2)[psi((r/10+1)/2) - psi(r/20)]

print("\nExact closed forms:")
n_sym = Symbol('n', integer=True, nonneg=True)

for r, label in [(1, "∫_0^∞"), (9, "∫_-∞^0")]:
    a = Rational(r, 10)
    # sum_{n>=0} (-1)^n / (10n+r) = (1/10) * sum_{n>=0} (-1)^n / (n + r/10)
    #   = (1/10) * (1/2)[psi((r/10+1)/2) - psi(r/20)]
    #   = (1/20)[psi((r+10)/20) - psi(r/20)]
    # so 5 * sum = (1/4)[psi((r+10)/20) - psi(r/20)]
    p_hi = Rational(r + 10, 20)
    p_lo = Rational(r, 20)
    exact = Rational(1, 4) * (polygamma(0, p_hi) - polygamma(0, p_lo))
    simplified = simplify(exact)
    print(f"  {label} = 5 * sum_n (-1)^n/(10n+{r})")
    print(f"         = (1/4)[psi({p_hi}) - psi({p_lo})]")
    print(f"         = {simplified}")
    print(f"         ≈ {float(simplified):.10f}")

print(f"\n  Total  = pi / (2*sin(pi/10)) ≈ {float(pi / (2 * sin(pi / 10))):.10f}")
