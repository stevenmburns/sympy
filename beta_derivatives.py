"""Compute high-order partial derivatives of the Beta function B(x, y)."""

from sympy import symbols, beta, diff, simplify, pprint

x, y = symbols('x y')

B = beta(x, y)

def beta_partial(order_x, order_y):
    """Return the partial derivative d^(m+n) B / dx^m dy^n."""
    return simplify(diff(B, (x, order_x), (y, order_y)))


if __name__ == "__main__":
    # Show derivatives up to 3rd order
    for m, n in ((2, 2), (3, 1)):
        expr = beta_partial(m, n)
        print(f"d^{m+n}B / dx^{m} dy^{n} =")
        pprint(expr)
        print()

    # Evaluate a specific high-order derivative at a point
    d = beta_partial(3, 1)
    print("d^4 B / dx^3 dy^1 evaluated at (1/4, 1/4):")
    print(float(d.subs({x: 0.25, y: 0.25})))
