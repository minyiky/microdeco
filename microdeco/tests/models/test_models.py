"""
Tests for the ZH_L16C_GF model in microdeco.models module.
"""

import pytest

from microdeco.test_utils import approx
from microdeco.models import ZH_L16C_GF, GradientFactors
from microdeco.gas import Gas
from microdeco import const


# --- Setup ---


def test_setup_compartments_count():
    """ZH-L16C should have exactly 16 compartments."""
    model = ZH_L16C_GF()
    assert len(model.compartments) == 16


def test_compartment_half_lives_order():
    """Nitrogen half-lives should be in ascending order (fast to slow tissues)."""
    model = ZH_L16C_GF()
    half_lives = [t.nitrogen_data.HALF_LIFE for t in model.compartments]
    assert half_lives == sorted(half_lives)
    assert half_lives[0] == 4.0
    assert half_lives[-1] == 635.0


def test_compartments_initial_state():
    """All compartments should start at surface air equilibrium."""
    model = ZH_L16C_GF()
    expected_pp_n2 = (1.0 - const.WATER_VAPOUR_PRESSURE) * const.PERCENTAGE_NITROGEN_IN_AIR

    for tissue in model.compartments:
        assert approx(tissue.pp_nitrogen, expected_pp_n2, rel=1e-9)
        assert tissue.pp_helium == 0.0


# --- load_tissues (constant depth) ---


def test_load_tissues_increases_nitrogen():
    """Loading at depth should increase nitrogen in all compartments."""
    model = ZH_L16C_GF()
    gas = Gas(p_o2=21, p_n2=79, p_he=0)

    initial_pp = [t.pp_nitrogen for t in model.compartments]
    model.load_tissues(30.0, 10.0, gas)

    for i, tissue in enumerate(model.compartments):
        assert tissue.pp_nitrogen > initial_pp[i]


def test_load_tissues_fast_compartment_loads_faster():
    """Compartment 1 (fast) should load more than compartment 16 (slow)."""
    model = ZH_L16C_GF()
    gas = Gas(p_o2=21, p_n2=79, p_he=0)
    initial_pp = model.compartments[0].pp_nitrogen

    model.load_tissues(30.0, 5.0, gas)

    fast_delta = model.compartments[0].pp_nitrogen - initial_pp
    slow_delta = model.compartments[-1].pp_nitrogen - initial_pp
    assert fast_delta > slow_delta


# --- load_tissues_with_rate ---


def test_load_tissues_with_rate_descent():
    """Tissues should load during descent."""
    model = ZH_L16C_GF()
    gas = Gas(p_o2=21, p_n2=79, p_he=0)

    initial_pp = [t.pp_nitrogen for t in model.compartments]
    model.load_tissues_with_rate(0.0, 30.0, 3.0, gas)

    for i, tissue in enumerate(model.compartments):
        assert tissue.pp_nitrogen > initial_pp[i]


def test_load_tissues_with_rate_helium():
    """Trimix descent should load both nitrogen and helium."""
    model = ZH_L16C_GF()
    trimix = Gas(p_o2=21, p_n2=35, p_he=44)

    model.load_tissues_with_rate(0.0, 40.0, 4.0, trimix)

    for tissue in model.compartments:
        assert tissue.pp_helium > 0.0


# --- calculate_ceilings ---


def test_ceiling_zero_at_surface():
    """Ceiling should be 0 (or negative, clamped to 0) at surface equilibrium."""
    model = ZH_L16C_GF()
    gf = GradientFactors(gf_low=0.3, gf_high=0.85)

    ceiling = model.calculate_ceilings(gf)
    assert ceiling == 0.0


def test_ceiling_positive_after_deep_dive():
    """After a deep dive, ceiling should be positive (diver needs deco)."""
    model = ZH_L16C_GF()
    gas = Gas(p_o2=21, p_n2=79, p_he=0)
    gf = GradientFactors(gf_low=0.3, gf_high=0.85)

    # Descend to 40m in 4 min, stay 30 min
    model.load_tissues_with_rate(0.0, 40.0, 4.0, gas)
    model.load_tissues(40.0, 30.0, gas)

    ceiling = model.calculate_ceilings(gf)
    assert ceiling > 0.0


def test_ceiling_returns_meters():
    """Ceiling should be in meters of seawater."""
    model = ZH_L16C_GF()
    gas = Gas(p_o2=21, p_n2=79, p_he=0)
    gf = GradientFactors(gf_low=0.3, gf_high=0.85)

    model.load_tissues_with_rate(0.0, 30.0, 3.0, gas)
    model.load_tissues(30.0, 20.0, gas)

    ceiling = model.calculate_ceilings(gf)
    # Ceiling for 30m/20min on air with GF 30/85 should be reasonable (0-20m range)
    assert 0.0 <= ceiling <= 20.0


def test_ceiling_higher_with_lower_gf():
    """More conservative GF (lower values) should produce a higher ceiling."""
    model = ZH_L16C_GF()
    gas = Gas(p_o2=21, p_n2=79, p_he=0)

    model.load_tissues_with_rate(0.0, 30.0, 3.0, gas)
    model.load_tissues(30.0, 20.0, gas)

    gf_conservative = GradientFactors(gf_low=0.2, gf_high=0.7)
    gf_liberal = GradientFactors(gf_low=0.5, gf_high=0.95)

    ceiling_conservative = model.calculate_ceilings(gf_conservative)
    ceiling_liberal = model.calculate_ceilings(gf_liberal)

    assert ceiling_conservative > ceiling_liberal


def test_ceiling_no_deco_shallow_dive():
    """A short shallow dive should not require deco (ceiling = 0)."""
    model = ZH_L16C_GF()
    gas = Gas(p_o2=21, p_n2=79, p_he=0)
    gf = GradientFactors(gf_low=0.3, gf_high=0.85)

    # 10m for 10 min - well within NDL even with conservative GF
    model.load_tissues_with_rate(0.0, 10.0, 1.0, gas)
    model.load_tissues(10.0, 10.0, gas)

    ceiling = model.calculate_ceilings(gf)
    assert ceiling == 0.0
