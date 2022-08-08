################################
# UNIT TESTING
################################
# Tests the functions in the Python module
#
# You can run the unit tests with [pytest](pytest.org)
# `pytest aav_ode_tests.py`

from re import T
import numpy as np

from aav_ode import Experiments

import pytest

import string


@pytest.fixture
def aav_model():

    aav_model = Experiments()
    return aav_model


def generate_random_string(length=20):  # pragma: no cover
    """Returns a random string of length"""

    return "".join([str(elem) for elem in np.random.choice([string.printable], length)])


def test_all():
    """Test end-to-end"""

    aav_model = Experiments()
    print(aav_model)
    t, prediction = aav_model.run_simulation()
    aav_model.plot_outputs(t, prediction)
    aav_model.plot_production(t, prediction)
    aav_model.plot_replication(t, prediction)

    np.testing.assert_almost_equal(np.max(t), 60.0)


def test_dna_fraction(aav_model):

    aav_model.run_simulation()
    np.testing.assert_almost_equal(aav_model.get_dna_fraction(), 32.99942077604155)


@pytest.mark.parametrize("k", [1e-5, 0.1, 1, 2, 3, 5, 10, 20, 25, 30])
def test_k(aav_model, k):

    for idx in range(len(aav_model.get_input_labels())):
        aav_model.set_k(idx, k)
        assert aav_model.get_k(idx) == k


def test_times(aav_model):

    t = 5
    aav_model.set_media_exchange_time(t)  # Hours
    aav_model.set_total_time(t)  # Hours

    assert aav_model.get_total_time() == t
    assert aav_model.get_media_exchange_time() == t


def test_calc_growth(aav_model):

    t, prediction = aav_model.run_simulation()

    np.testing.assert_almost_equal(
        np.max(aav_model._calc_growth_rate(t)), 0.0494279581111241
    )


def test_convert_ml(aav_model):

    t, prediction = aav_model.run_simulation()

    np.testing.assert_almost_equal(
        aav_model._convert_to_ml(t, prediction)[1, 1], 48879421341.434494,
        decimal=3
    )

def test_get_input_labels(aav_model):

    l = aav_model.get_input_labels()

    assert len(l) == 14
    assert l[0] == "kUptake"


def test_get_initial_concentrations(aav_model):

    l = aav_model.get_initial_concentrations()

    assert len(l) == 3
    assert l[0] == 76000.0


def test_get_solver_descriptions(aav_model):

    d = aav_model.get_solver_description()

    assert d["RK45"] == "Explicit Runge-Kutta method of order 5"
    assert (
        d["BDF"] == "Implicit multi-step variable-order (1 to 5) "
        "method based on a backward differentiation formula "
        "for the derivative approximation "
    )


def test_get_set_solver(aav_model):

    aav_model.set_ode_solver_method("BDF")
    assert aav_model.get_ode_solver_method() == "BDF"


def test_get_set_output_label(aav_model):

    l = generate_random_string(length=20)
    aav_model.set_output_label(pRC_extracellular=l)

    assert aav_model.get_output_label(0) == l

    assert aav_model.get_output_labels()[0] == l
