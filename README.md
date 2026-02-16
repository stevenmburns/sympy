# Beta Function Partial Derivatives

Symbolic computation of high-order partial derivatives of the [Beta function](https://en.wikipedia.org/wiki/Beta_function) B(x, y) = Γ(x)Γ(y) / Γ(x+y) using SymPy.

## Setup

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install sympy
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

## Output

Derivatives are expressed symbolically in terms of polygamma functions and B(x, y), e.g.:

```
d^4B / dx^3 dy^1 evaluated at (1/4, 1/4):
2.8773310532591805
```
