# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Symbolic computation of high-order partial derivatives of the Beta function B(x, y) using SymPy.

## Setup

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Running

```bash
.venv/bin/python beta_derivatives.py
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

Single-file project (`beta_derivatives.py`):
- Module-level symbols `x`, `y` and `B = beta(x, y)` are shared across the module
- `beta_partial(order_x, order_y)` — computes `∂^(m+n) B / ∂x^m ∂y^n` using `sympy.diff` and `simplify`
- Results are expressed in terms of polygamma functions times B(x, y)
- Can be imported as a library or run as a script
