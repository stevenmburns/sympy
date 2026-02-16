"""Tests for beta_partial derivatives."""

from beta_derivatives import beta_partial, x, y


def test_d4_dx2_dy2():
    d = beta_partial(2, 2)
    result = float(d.subs({x: 0.25, y: 0.25}))
    assert abs(result - 0.4965) < 0.001


def test_d4_dx3_dy1():
    d = beta_partial(3, 1)
    result = float(d.subs({x: 0.25, y: 0.25}))
    assert abs(result - 2.877) < 0.001
