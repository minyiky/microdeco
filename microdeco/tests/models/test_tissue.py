"""
Tests for the Tissue and GasCoeffients classes in microdeco.models module.
"""

import pytest

from microdeco.test_utils import approx
from microdeco.models import Tissue, GasCoeffients, GradientFactors
from microdeco import const


# --- GasCoeffients tests ---


def test_gas_coeffients_k_calculation():
    """Verify k is correctly derived from HALF_LIFE."""
    gc = GasCoeffients(A=1.2599, B=0.5050, HALF_LIFE=4.0)
    assert approx(gc.k, const.LOG_2 / 4.0, rel=1e-9)


def test_gas_coeffients_stores_values():
    """Verify A, B, HALF_LIFE are stored correctly."""
    gc = GasCoeffients(A=0.5043, B=0.8434, HALF_LIFE=38.3)
    assert gc.A == 0.5043
    assert gc.B == 0.8434
    assert gc.HALF_LIFE == 38.3


# --- Tissue init tests ---


def test_tissue_initial_pp_nitrogen():
    """Verify initial nitrogen partial pressure assumes air at surface."""
    n2 = GasCoeffients(A=1.2599, B=0.5050, HALF_LIFE=4.0)
    he = GasCoeffients(A=1.7424, B=0.4245, HALF_LIFE=1.51)
    tissue = Tissue(nitrogen=n2, helium=he)

    expected_pp_n2 = (1.0 - const.WATER_VAPOUR_PRESSURE) * const.PERCENTAGE_NITROGEN_IN_AIR
    assert approx(tissue.pp_nitrogen, expected_pp_n2, rel=1e-9)
    assert tissue.pp_helium == 0.0


def test_tissue_initial_pp_with_custom_pressure():
    """Verify initial pp_nitrogen scales with initial_pressure."""
    n2 = GasCoeffients(A=1.2599, B=0.5050, HALF_LIFE=4.0)
    he = GasCoeffients(A=1.7424, B=0.4245, HALF_LIFE=1.51)
    tissue = Tissue(nitrogen=n2, helium=he, initial_pressure=2.0)

    expected_pp_n2 = (2.0 - const.WATER_VAPOUR_PRESSURE) * const.PERCENTAGE_NITROGEN_IN_AIR
    assert approx(tissue.pp_nitrogen, expected_pp_n2, rel=1e-9)


# --- Tissue load tests (constant depth, rate=0) ---


def test_tissue_load_constant_depth_nitrogen_increases():
    """Loading at depth should increase nitrogen partial pressure."""
    n2 = GasCoeffients(A=1.2599, B=0.5050, HALF_LIFE=4.0)
    he = GasCoeffients(A=1.7424, B=0.4245, HALF_LIFE=1.51)
    tissue = Tissue(nitrogen=n2, helium=he)

    initial_pp = tissue.pp_nitrogen
    # Simulate 30m depth on air: abs_p = 4.0, p_n2 = 0.79 * (4.0 - 0.0627)
    p_n2 = 0.79 * (4.0 - const.WATER_VAPOUR_PRESSURE)
    tissue.load(p_n2=p_n2, p_he=0.0, t=10.0)

    assert tissue.pp_nitrogen > initial_pp


def test_tissue_load_constant_depth_nitrogen_converges():
    """After long time at depth, tissue pressure should approach alveolar pressure."""
    n2 = GasCoeffients(A=1.2599, B=0.5050, HALF_LIFE=4.0)
    he = GasCoeffients(A=1.7424, B=0.4245, HALF_LIFE=1.51)
    tissue = Tissue(nitrogen=n2, helium=he)

    # 30m depth, air
    p_n2 = 0.79 * (4.0 - const.WATER_VAPOUR_PRESSURE)
    # 40 minutes = 10 half-lives for this compartment, should be ~99.9% saturated
    tissue.load(p_n2=p_n2, p_he=0.0, t=40.0)

    assert approx(tissue.pp_nitrogen, p_n2, rel=1e-2)


# --- Tissue load tests (variable depth, rate!=0) ---


def test_tissue_load_with_rate_descent():
    """Loading during descent should increase nitrogen partial pressure."""
    n2 = GasCoeffients(A=1.2599, B=0.5050, HALF_LIFE=4.0)
    he = GasCoeffients(A=1.7424, B=0.4245, HALF_LIFE=1.51)
    tissue = Tissue(nitrogen=n2, helium=he)

    initial_pp = tissue.pp_nitrogen
    # Descending from 0 to 30m in 3 min: rate = (4.0-1.0)/3 = 1.0 bar/min
    p_n2 = 0.79 * (1.0 - const.WATER_VAPOUR_PRESSURE)  # starting depth alveolar pressure
    tissue.load(p_n2=p_n2, p_he=0.0, t=3.0, rate=0.79 * 1.0)

    assert tissue.pp_nitrogen > initial_pp


def test_tissue_load_with_rate_ascent():
    """Loading during ascent should decrease nitrogen partial pressure (if saturated)."""
    n2 = GasCoeffients(A=1.2599, B=0.5050, HALF_LIFE=4.0)
    he = GasCoeffients(A=1.7424, B=0.4245, HALF_LIFE=1.51)
    tissue = Tissue(nitrogen=n2, helium=he)

    # First saturate at 30m
    p_n2_30m = 0.79 * (4.0 - const.WATER_VAPOUR_PRESSURE)
    tissue.load(p_n2=p_n2_30m, p_he=0.0, t=40.0)
    pp_at_depth = tissue.pp_nitrogen

    # Now ascend from 30m to 0m in 3 min
    p_n2_surface = 0.79 * (4.0 - const.WATER_VAPOUR_PRESSURE)  # starting from 30m
    rate = 0.79 * ((1.0 - 4.0) / 3.0)  # negative rate for ascent
    tissue.load(p_n2=p_n2_surface, p_he=0.0, t=3.0, rate=rate)

    assert tissue.pp_nitrogen < pp_at_depth


# --- Tissue helium loading tests ---


def test_tissue_load_helium_from_zero():
    """Helium should start loading even when pp_helium is initially 0."""
    n2 = GasCoeffients(A=1.2599, B=0.5050, HALF_LIFE=4.0)
    he = GasCoeffients(A=1.7424, B=0.4245, HALF_LIFE=1.51)
    tissue = Tissue(nitrogen=n2, helium=he)

    assert tissue.pp_helium == 0.0

    # Breathe trimix 21/35 (44% He) at 30m
    abs_p = 4.0
    p_he = 0.44 * (abs_p - const.WATER_VAPOUR_PRESSURE)
    p_n2 = 0.35 * (abs_p - const.WATER_VAPOUR_PRESSURE)
    tissue.load(p_n2=p_n2, p_he=p_he, t=10.0)

    assert tissue.pp_helium > 0.0


def test_tissue_load_helium_with_rate():
    """Helium should load during descent on trimix."""
    n2 = GasCoeffients(A=1.2599, B=0.5050, HALF_LIFE=4.0)
    he = GasCoeffients(A=1.7424, B=0.4245, HALF_LIFE=1.51)
    tissue = Tissue(nitrogen=n2, helium=he)

    # Descend to 30m on trimix 21/35
    p_he = 0.44 * (1.0 - const.WATER_VAPOUR_PRESSURE)
    p_n2 = 0.35 * (1.0 - const.WATER_VAPOUR_PRESSURE)
    rate = (4.0 - 1.0) / 3.0  # 1 bar/min
    tissue.load(p_n2=p_n2, p_he=p_he, t=3.0, rate=rate)

    assert tissue.pp_helium > 0.0


# --- Tissue compute_ceiling tests ---


def test_tissue_ceiling_at_surface():
    """A tissue at surface equilibrium should have a negative or zero ceiling."""
    n2 = GasCoeffients(A=1.2599, B=0.5050, HALF_LIFE=4.0)
    he = GasCoeffients(A=1.7424, B=0.4245, HALF_LIFE=1.51)
    tissue = Tissue(nitrogen=n2, helium=he)
    gf = GradientFactors(gf_low=0.3, gf_high=0.85)

    ceiling = tissue.compute_ceiling(gf)
    # Ceiling in bar - at surface equilibrium it should be below surface pressure
    assert ceiling < const.SURFACE_PRESSURE


def test_tissue_ceiling_increases_with_depth_exposure():
    """Ceiling should increase after exposure at depth."""
    n2 = GasCoeffients(A=1.2599, B=0.5050, HALF_LIFE=4.0)
    he = GasCoeffients(A=1.7424, B=0.4245, HALF_LIFE=1.51)
    tissue = Tissue(nitrogen=n2, helium=he)
    gf = GradientFactors(gf_low=0.3, gf_high=0.85)

    ceiling_before = tissue.compute_ceiling(gf)

    # Saturate at 30m on air
    p_n2 = 0.79 * (4.0 - const.WATER_VAPOUR_PRESSURE)
    tissue.load(p_n2=p_n2, p_he=0.0, t=20.0)

    ceiling_after = tissue.compute_ceiling(gf)
    assert ceiling_after > ceiling_before


# --- GradientFactors valid cases ---


def test_gradient_factors_defaults():
    """Verify default GF values."""
    gf = GradientFactors()
    assert gf.gf_low == 0.3
    assert gf.gf_high == 0.85
    assert gf.initial_depth == 0.0


def test_gradient_factors_set_initial_depth():
    """Verify set_initial_depth stores the value."""
    gf = GradientFactors(gf_low=0.3, gf_high=0.85)
    gf.set_initial_depth(30.0)
    assert gf.initial_depth == 30.0


def test_gradient_factors_get_current_gf_at_depth():
    """At depths deeper than initial, GF should be gf_low."""
    gf = GradientFactors(gf_low=0.3, gf_high=0.85)
    gf.set_initial_depth(20.0)
    assert gf.get_current_gf(25.0) == 0.3


def test_gradient_factors_get_current_gf_at_surface():
    """At surface (0m), GF should be gf_high."""
    gf = GradientFactors(gf_low=0.3, gf_high=0.85)
    gf.set_initial_depth(20.0)
    assert approx(gf.get_current_gf(0.0), 0.85, rel=1e-9)


def test_gradient_factors_get_current_gf_interpolation():
    """At half the initial depth, GF should be midpoint of gf_low and gf_high."""
    gf = GradientFactors(gf_low=0.3, gf_high=0.85)
    gf.set_initial_depth(20.0)
    expected = 0.3 + (0.85 - 0.3) * 0.5  # 0.575
    assert approx(gf.get_current_gf(10.0), expected, rel=1e-9)
