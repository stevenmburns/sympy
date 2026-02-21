# Beta Function Partial Derivatives

[![CI](https://github.com/stevenmburns/sympy/actions/workflows/ci.yml/badge.svg)](https://github.com/stevenmburns/sympy/actions/workflows/ci.yml)

Symbolic computation of high-order partial derivatives of the [Beta function](https://en.wikipedia.org/wiki/Beta_function) B(x, y) = Γ(x)Γ(y) / Γ(x+y) using SymPy.

## Setup

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Usage

Run the script directly to see example derivatives:

```bash
python beta_derivatives.py
```

Or import `beta_partial` for custom use:

```python
from beta_derivatives import beta_partial, x, y

# Compute ∂⁵B / ∂x³∂y²
expr = beta_partial(3, 2)

# Evaluate numerically
result = float(expr.subs({x: 2, y: 3}))
```

## Testing

```bash
pytest test_beta_derivatives.py -v
```

## Output

Derivatives are expressed symbolically in terms of polygamma functions and B(x, y), e.g.:

```
d^4B / dx^3 dy^1 evaluated at (1/4, 1/4):
2.8773310532591805
```

---

## Numerical Integration: `numerical_integration.py`

Numerical integration of `exp(4x/5 − ln 2) · sech(x)` over the real line using SciPy.

### The integral

```
∫_{-∞}^{∞} exp(4x/5 − ln 2) · sech(x) dx  =  π / (2 sin(π/10))  ≈  5.08320369...
```

### Numerically stable evaluation

The integrand `f(x) = exp(4x/5) / (eˣ + e⁻ˣ)` is rewritten to avoid overflow:

- x ≥ 0: `f(x) = exp(−x/5) / (1 + exp(−2x))`
- x < 0: `f(x) = exp(9x/5) / (1 + exp(2x))`

### Asymptotic expansions

Expanding `1/(1 + small)` as a geometric series in each half:

**x → +∞:**
```
f(x) = Σ_{n≥0} (−1)ⁿ exp(−(2n + 1/5)x)
     = exp(−x/5) − exp(−11x/5) + exp(−21x/5) − ···
```

**x → −∞:**
```
f(x) = Σ_{n≥0} (−1)ⁿ exp((2n + 9/5)x)
     = exp(9x/5) − exp(19x/5) + exp(29x/5) − ···
```

These are not merely asymptotic: each series converges for all x > 0 (resp. x < 0).

### Integrating the series term by term

```
∫_0^∞  f(x) dx  =  Σ_{n≥0} (−1)ⁿ / (2n + 1/5)  =  5 Σ_{n≥0} (−1)ⁿ / (10n + 1)
∫_{-∞}^0 f(x) dx  =  Σ_{n≥0} (−1)ⁿ / (2n + 9/5)  =  5 Σ_{n≥0} (−1)ⁿ / (10n + 9)
```

### Exact closed forms via the digamma function

Using `Σ_{n≥0} (−1)ⁿ / (n + a) = (1/2)[ψ((a+1)/2) − ψ(a/2)]`:

```
∫_0^∞     = (1/4)[ψ(11/20) − ψ( 1/20)]  ≈  4.69047143516...
∫_{-∞}^0  = (1/4)[ψ(19/20) − ψ( 9/20)]  ≈  0.39273225715...
─────────────────────────────────────────────────────────────
Total      = π / (2 sin(π/10))            ≈  5.08320369231...
```

Verified to 30 significant figures: the two digamma differences sum exactly to `π / (2 sin(π/10))`.
