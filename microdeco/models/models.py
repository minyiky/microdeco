from .tissue import GasCoeffients, Tissue
from ..gas import Gas
from .gradient_factor import GradientFactors


class ZH_GF:
    pass


class ZH_L16C_GF:  # source: ostc firmware code
    """
    ZH-L16C-GF decompression model.
    """

    @staticmethod
    def setup_compartments() -> list[Tissue]:
        return [
            Tissue(
                nitrogen=GasCoeffients(A=1.2599, B=0.5050, HALF_LIFE=4.0),
                helium=GasCoeffients(A=1.7424, B=0.4245, HALF_LIFE=1.51),
            ),
            Tissue(
                nitrogen=GasCoeffients(A=1.0000, B=0.6514, HALF_LIFE=8.0),
                helium=GasCoeffients(A=1.3830, B=0.5747, HALF_LIFE=3.02),
            ),
            Tissue(
                nitrogen=GasCoeffients(A=0.8618, B=0.7222, HALF_LIFE=12.5),
                helium=GasCoeffients(A=1.1919, B=0.6527, HALF_LIFE=4.72),
            ),
            Tissue(
                nitrogen=GasCoeffients(A=0.7562, B=0.7825, HALF_LIFE=18.5),
                helium=GasCoeffients(A=1.0458, B=0.7223, HALF_LIFE=6.99),
            ),
            Tissue(
                nitrogen=GasCoeffients(A=0.6200, B=0.8126, HALF_LIFE=27.0),
                helium=GasCoeffients(A=0.9220, B=0.7582, HALF_LIFE=10.21),
            ),
            Tissue(
                nitrogen=GasCoeffients(A=0.5043, B=0.8434, HALF_LIFE=38.3),
                helium=GasCoeffients(A=0.8205, B=0.7957, HALF_LIFE=14.48),
            ),
            Tissue(
                nitrogen=GasCoeffients(A=0.4410, B=0.8693, HALF_LIFE=54.3),
                helium=GasCoeffients(A=0.7305, B=0.8279, HALF_LIFE=20.53),
            ),
            Tissue(
                nitrogen=GasCoeffients(A=0.4000, B=0.8910, HALF_LIFE=77.0),
                helium=GasCoeffients(A=0.6502, B=0.8553, HALF_LIFE=29.11),
            ),
            Tissue(
                nitrogen=GasCoeffients(A=0.3750, B=0.9092, HALF_LIFE=109.0),
                helium=GasCoeffients(A=0.5950, B=0.8757, HALF_LIFE=41.20),
            ),
            Tissue(
                nitrogen=GasCoeffients(A=0.3500, B=0.9222, HALF_LIFE=146.0),
                helium=GasCoeffients(A=0.5545, B=0.8903, HALF_LIFE=55.19),
            ),
            Tissue(
                nitrogen=GasCoeffients(A=0.3295, B=0.9319, HALF_LIFE=187.0),
                helium=GasCoeffients(A=0.5333, B=0.8997, HALF_LIFE=70.69),
            ),
            Tissue(
                nitrogen=GasCoeffients(A=0.3065, B=0.9403, HALF_LIFE=239.0),
                helium=GasCoeffients(A=0.5189, B=0.9073, HALF_LIFE=90.34),
            ),
            Tissue(
                nitrogen=GasCoeffients(A=0.2835, B=0.9477, HALF_LIFE=305.0),
                helium=GasCoeffients(A=0.5181, B=0.9122, HALF_LIFE=115.29),
            ),
            Tissue(
                nitrogen=GasCoeffients(A=0.2610, B=0.9544, HALF_LIFE=390.0),
                helium=GasCoeffients(A=0.5176, B=0.9171, HALF_LIFE=147.42),
            ),
            Tissue(
                nitrogen=GasCoeffients(A=0.2480, B=0.9602, HALF_LIFE=498.0),
                helium=GasCoeffients(A=0.5172, B=0.9217, HALF_LIFE=188.24),
            ),
            Tissue(
                nitrogen=GasCoeffients(A=0.2327, B=0.9653, HALF_LIFE=635.0),
                helium=GasCoeffients(A=0.5119, B=0.9267, HALF_LIFE=240.03),
            ),
        ]

    def __init__(self):
        self.compartments = self.setup_compartments()

    def load_tissues(self, depth: float, t: float, gas: Gas) -> None:
        """
        Load all tissues in the model with the given depth, time, and gas mixture.

        This is a convenience method for constant depth segments (rate = 0).
        For ascending or descending profiles, use load_tissues_with_rate().

        Args:
            depth (float): The current depth in meters.
            t (float): The time in minutes.
            gas (Gas): The breathing gas mixture.
        """
        from .. import const
        abs_p = depth / 10 + 1
        p_n2 = gas.p_n2 * (abs_p - const.WATER_VAPOUR_PRESSURE)
        p_he = gas.p_he * (abs_p - const.WATER_VAPOUR_PRESSURE)
        for tissue in self.compartments:
            tissue.load(p_n2=p_n2, p_he=p_he, t=t)

    def load_tissues_with_rate(self, start_depth: float, end_depth: float, t: float, gas: Gas) -> None:
        """
        Load all tissues in the model while changing depth at a constant rate.

        This uses the Schreiner equation to account for the changing ambient pressure
        during ascent or descent. The alveolar partial pressure is calculated from
        the starting depth, and the pressure rate is derived from the depth change.

        Args:
            start_depth (float): Starting depth in meters.
            end_depth (float): Ending depth in meters.
            t (float): The time in minutes to transition from start to end depth.
            gas (Gas): The breathing gas mixture.
        """
        from .. import const

        # Calculate absolute pressures
        start_abs_p = start_depth / 10 + 1
        end_abs_p = end_depth / 10 + 1

        # Pressure change per minute in bar/min (positive for descent, negative for ascent)
        rate = (end_abs_p - start_abs_p) / t

        # Alveolar partial pressure calculated from starting depth (following decotengu approach)
        p_n2 = gas.p_n2 * (start_abs_p - const.WATER_VAPOUR_PRESSURE)
        p_he = gas.p_he * (start_abs_p - const.WATER_VAPOUR_PRESSURE)

        for tissue in self.compartments:
            tissue.load(p_n2=p_n2, p_he=p_he, t=t, rate=rate)

    def calculate_ceilings(self, gf: GradientFactors) -> float:
        """
        Calculate the ceiling for all tissues in the model using the given gradient factors.
        Args:
            gf (GradientFactors): The gradient factors to use for the calculation.
        Returns:
            float: The maximum ceiling pressure of all tissues in meters.
        """
        ceiling = max([tissue.compute_ceiling(gf) for tissue in self.compartments])
        return max((ceiling - 1) * 10, 0)  # convert bar to meters of sea water
