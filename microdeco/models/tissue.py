from .. import const
from ..equations import load_tissue_constant, load_tissue, buhlmann
from . import GradientFactors


class Tissue(object):
    """
    A class to represent a tissue.

    Attributes:
        nitrogen (GasCoeffients): The gas coefficients for Nitrogen.
        helium (GasCoeffients): The gas coefficients for Helium.
    """
    __slots__ = ('nitrogen_data', 'helium_data', 'pp_nitrogen', 'pp_helium')

    def __init__(
        self,
        nitrogen: "GasCoeffients",
        helium: "GasCoeffients",
        initial_pressure: float = 1.0,
    ):
        """
        Initialize a new Tissue object.

        Args:
            nitrogen (GasCoeffients): The gas coefficients for Nitrogen.
            helium (GasCoeffients): The gas coefficients for Helium.
        """
        self.nitrogen_data = nitrogen
        self.helium_data = helium

        # Initial partial pressures in the tissue using air as the starting gas
        self.pp_nitrogen = (
            (initial_pressure - const.WATER_VAPOUR_PRESSURE) * const.PERCENTAGE_NITROGEN_IN_AIR
        )  # assuming air
        self.pp_helium = 0.0

    def load(self, p_n2: float, p_he: float, t: float, rate: float = 0.0) -> None:
        """
        Load the tissue with the given partial pressures of Nitrogen and Helium over time.

        This function updates the partial pressures of Nitrogen and Helium in the tissue
        using the Schreiner equation (for non-zero rate) or Haldane equation (for zero rate).

        Args:
            p_n2 (float): The partial pressure of Nitrogen.
            p_he (float): The partial pressure of Helium.
            t (float): The time in minutes.
            rate (float): The rate of change of ambient pressure in bar/min.
                         Positive for descent, negative for ascent, 0 for constant depth.
        """
        if rate == 0:
            self.pp_nitrogen = load_tissue_constant(
                p_alv=p_n2,
                p_i=self.pp_nitrogen,
                k=self.nitrogen_data.k,
                t=t,
            )
            if p_he > 0.0 or self.pp_helium > 0.0:
                self.pp_helium = load_tissue_constant(
                    p_alv=p_he,
                    p_i=self.pp_helium,
                    k=self.helium_data.k,
                    t=t,
                )
        else:
            self.pp_nitrogen = load_tissue(
                p_alv=p_n2,
                p_i=self.pp_nitrogen,
                k=self.nitrogen_data.k,
                t=t,
                R=rate,
            )
            if p_he > 0.0 or self.pp_helium > 0.0:
                self.pp_helium = load_tissue(
                    p_alv=p_he,
                    p_i=self.pp_helium,
                    k=self.helium_data.k,
                    t=t,
                    R=rate,
                )

    def compute_ceiling(self, gf) -> float:
        """
        Compute the ceiling pressure for the tissue.

        Accepts either a float GF value directly (used during deco schedule
        calculation with interpolated GFs) or a GradientFactors object (legacy,
        uses gf_low).

        Args:
            gf: A float gradient factor (0-1) or a GradientFactors object.
        Returns:
            float: The ceiling pressure for the tissue in bar.
        """
        gf_value = gf if isinstance(gf, float) else gf.gf_low
        return buhlmann(
            gf=gf_value,
            p_n2=self.pp_nitrogen,
            p_he=self.pp_helium,
            A_n2=self.nitrogen_data.A,
            A_he=self.helium_data.A,
            B_n2=self.nitrogen_data.B,
            B_he=self.helium_data.B,
        )


class GasCoeffients:
    """
    A class to represent a set of gas coefficients for a tissue.

    Attributes:
        A (float): The A coefficient for Nitrogen.
        B (float): The B coefficient for Nitrogen.
        HALF_LIFE (float): The half-life of the tissue for Nitrogen.
        k (float): The gas decay constant for the tissue.
    """
    __slots__ = ('A', 'B', 'HALF_LIFE', 'k')

    def __init__(self, A: float, B: float, HALF_LIFE: float):
        """
        Initialize a new GasCoeffients object.

        Args:
            A (float): The A coefficient for Nitrogen.
            B (float): The B coefficient for Nitrogen.
            HALF_LIFE (float): The half-life of the tissue for Nitrogen.
        """
        self.A = A
        self.B = B
        self.HALF_LIFE = HALF_LIFE
        self.k = const.LOG_2 / HALF_LIFE
