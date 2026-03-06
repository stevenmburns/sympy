# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Symbolic and numerical computation using SymPy, NumPy, and SciPy. Includes:
- High-order partial derivatives of the Beta function B(x, y)
- Numerical integration of definite integrals with series derivations and convergence acceleration

## Setup

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Running

```bash
.venv/bin/python beta_derivatives.py
.venv/bin/python numerical_integration.py
.venv/bin/python arctan_log_integral.py
.venv/bin/python erf_integral.py
.venv/bin/python plot.py
```

## Testing

```bash
.venv/bin/python -m pytest test_beta_derivatives.py -v
```

## Linting

```bash
ruff check .
```

## Architecture

**`beta_derivatives.py`** — Beta function partial derivatives:
- Module-level symbols `x`, `y` and `B = beta(x, y)` are shared across the module
- `beta_partial(order_x, order_y)` — computes `∂^(m+n) B / ∂x^m ∂y^n` using `sympy.diff` and `simplify`
- Results are expressed in terms of polygamma functions times B(x, y)
- Can be imported as a library or run as a script

**`numerical_integration.py`** — Numerical integration of `exp(4x/5 - ln2)·sech(x)` over ℝ:
- Numerically stable integrand, asymptotic expansions, term-by-term series integration
- Exact closed form `π/(2·sin(π/10))` derived via digamma functions

**`arctan_log_integral.py`** — Numerical integration of `arctan(x)/x · ln(1+x)` over [0,1]:
- Two digamma series derived by expanding either factor as a power series and integrating term-by-term
- Each term reduces via IBP + substitution `t=x²` using `∫₀¹ t^{s-1}/(1+t) dt = ½[ψ((s+1)/2) - ψ(s/2)]`
- Series 1: `I = Σ_{k≥1} (-1)^{k-1}/(4k²) · [π + ψ((k+1)/4) - ψ((k+3)/4)]`
- Series 2: `I = Σ_{n≥0} (-1)^n/(2n+1)² · [ln2 - ½(ψ(n+3/2) - ψ(n+1))]`
- Iterated Aitken Δ² acceleration reaches machine precision in 4 passes from 50 raw terms
- Euler transform does NOT work here (terms not monotonically decreasing); use Aitken

**`erf_integral.py`** — Numerical integration of `e^{-x} sin(x) erf(x)` over [0,∞):
- IBP with `u = erf(x)` reduces to `(1/√π) ∫₀^∞ e^{-x²-x}(sin x + cos x) dx`
- Completing the square + complex Gaussian integral identity gives closed form `(1/√2) Im[e^{i(π/4-1/2)} erfc(1/2 - i/2)]`
- `scipy.special.erf` accepts complex arguments; `erfc(z) = 1 - erf(z)` used directly

**`plot.py`** — Plots `g(x) = ||x³-x| - x²|` over [0,2] with piecewise polynomial fit:
- Endpoints at golden ratio `φ = (1+√5)/2`; exact integrals expressed in terms of √5
- Requires `from sympy import sqrt, simplify, Rational`

## Script conventions for new integrals

Each integral script follows a three-layer structure:

1. **Direct numerical reference** — `scipy.integrate.quad` on the raw integrand first. This is the ground truth everything else is checked against.

2. **Analytic reduction** — reduce the integral to a simpler form by hand (IBP, series expansion, substitution, etc.) and verify numerically that the reduced form agrees with layer 1. Show the derivation in the module docstring.

3. **Closed form** — express the result in terms of named special functions (`polygamma`, `erfc` at complex argument, etc.) and verify it matches to machine precision. Print the discrepancy explicitly.

**Common reduction techniques used in this repo:**
- IBP to strip off a special-function factor (erf, arctan, ln) and expose a Gaussian or exponential kernel
- Expanding a factor as a power series and integrating term-by-term, then identifying each term via a known integral identity
- Completing the square to convert real exponentials to complex Gaussian form, then using `∫_a^∞ e^{-t²+ibt} dt = e^{-b²/4}·(√π/2) erfc(a - ib/2)`
- Digamma reduction: `∫₀¹ t^{s-1}/(1+t) dt = ½[ψ((s+1)/2) - ψ(s/2)]`

**Convergence acceleration:** for slowly converging alternating series, use iterated Aitken Δ² (not Euler transform — Euler requires monotonically decreasing terms).
