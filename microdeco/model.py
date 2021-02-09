# =============================================================================

# 

# =============================================================================



try:
    from ucollections import namedtuple
except ImportError:
    from collections import namedtuple
from math import exp



from . import const

from .algorithms import eq_gf_limit



Data = namedtuple('Data', 'tissues gf')



class ZH_L16_GF():

    """

    Base abstract class for Buhlmann ZH-L16 decompression model with

    gradient factors by Erik Baker - ZH-L16B-GF.

    """



    NUM_COMPARTMENTS = 16

    N2_A = None

    N2_B = None

    HE_A = None

    HE_B = None

    N2_HALF_LIFE = None

    HE_HALF_LIFE = None

    START_P_N2 = 0.7902 # starting pressure of N2 in tissues

    START_P_HE = 0.0    # starting pressure of He in tissues



    def __init__(self, gf_low=0.3, gf_high=0.8):

        """

        Create instance of the model.



        :var gf_low: Gradient factor low parameter.

        :var gf_high: Gradient factor high parameter.

        """

        self.n2_k_const = self._k_const(self.N2_HALF_LIFE)

        self.he_k_const = self._k_const(self.HE_HALF_LIFE)

        self.gf_low = gf_low

        self.gf_high = gf_high

        self.water_vapour_pressure = const.WATER_VAPOUR_PRESSURE_DEFAULT





    def init(self, surface_pressure):

        """

        Initialize pressure of inert gas in all tissues.



        The method uses starting tissue pressure values for nitrogen and

        helium.



        :param surface_pressure: Surface pressure [bar].

        """

        p_n2 = self.START_P_N2 * (surface_pressure - self.water_vapour_pressure)

        p_he = self.START_P_HE

        data = Data(tuple([(p_n2, p_he)] * self.NUM_COMPARTMENTS), self.gf_low)

        return data





    def load(self, abs_p, time, gas, rate, data):

        """

        Calculate gas loading for all tissue compartments.



        The method returns decompression data model information.



        :param abs_p: Absolute pressure [bar] (current depth).

        :param time: Time of exposure [min] (i.e. time of ascent).

        :param gas: Gas mix configuration.

        :param rate: Pressure rate change [bar/min].

        :param data: Decompression model data.



        .. seealso::



            - :py:meth:`decotengu.model.ZH_L16_GF._tissue_loaders`

            - :py:meth:`decotengu.model.ZH_L16_GF._tissue_loader`

        """

        n2_loader, he_loader = self._tissue_loaders(abs_p, gas, rate)



        tp = tuple(

            (n2_loader(time, p_n2, i), he_loader(time, p_he, i))

            for i, (p_n2, p_he) in enumerate(data.tissues)

        )

        return Data(tp, data.gf)





    def ceiling_limit(self, data, gf=None):

        """

        Calculate pressure of ascent ceiling limit using decompression

        model data.



        The pressure is the shallowest depth a diver can reach without

        decompression sickness. If pressure limit is 3 bar, then diver

        should not go shallower than 20m.



        FIXME: the method signature is gradient factor specific, the

            signature has to be made decompression model independent



        :param data: Decompression model data.

        :param gf: Gradient factor value, `gf_low` by default.



        .. seealso::



            - :py:meth:`decotengu.model.ZH_L16_GF.gf_limit`

            - :py:meth:`decotengu.model.ZH_L16_GF._tissue_loader`

        """

        return max(self.gf_limit(gf, data))





    def _k_const(self, half_life):

        """

        Calculate gas decay constant :math:`k` for each tissue compartment

        half-life value.



        :param half_life: Collection of half-life values for each tissue

            compartment.

        """

        return tuple(const.LOG_2 / v for v in half_life)





    def _tissue_loaders(self, abs_p, gas, rate):

        """

        Create function to load tissue compartment with inert gas for each

        inert gas specified in gas mix configuration.



        :param abs_p: Absolute pressure of current depth [bar] (:math:`P_{abs}`).

        :param gas: Gas mix configuration.

        :param rate: Pressure rate change [bar/min] (:math:`P_{rate}`).

        """

        n2_loader = self._tissue_loader(

            abs_p, gas.n2 / 100, rate, self.n2_k_const

        )

        he_loader = self._tissue_loader(

            abs_p, gas.he / 100, rate, self.he_k_const

        )

        return n2_loader, he_loader





    def _tissue_loader(self, abs_p, f_gas, rate, k_const):

        """

        Create function to load tissue compartment with inert gas.



        The created function uses Schreiner equation and has the following

        parameters



        time

            Time of exposure [min] at depth (:math:`T_{time}`).

        p_i

            Initial (current) pressure of inert gas in tissue compartment

            [bar] (:math:`P_{i}`).

        tissue_no

            Number of tissue compartment in the decompression model

            (starting with zero).



        See :ref:`eq-schreiner` section for details.



        :param abs_p: Absolute pressure of current depth [bar] (:math:`P_{abs}`).

        :param f_gas: Inert gas fraction, i.e. for air it is 0.79 (:math:`F_{gas}`).

        :param rate: Pressure rate change [bar/min] (:math:`P_{rate}`).

        :param k_const: Collection of gas decay constants for each tissue

            compartment (:math:`k`).

        """

        p_alv = f_gas * (abs_p - self.water_vapour_pressure)

        r = f_gas * rate

        def f(time, p_i, tissue_no):

            assert time > 0

            k = k_const[tissue_no]

            return p_alv + r * (time - 1 / k) - (p_alv - p_i - r / k) \
                * exp(-k * time)

        return f





    def gf_limit(self, gf, data):

        """

        Calculate pressure of ascent ceiling for each tissue compartment.



        The method returns a tuple of values - a pressure value for each

        tissue compartment.



        :param gf: Gradient factor.

        :param data: Decompression model data.

        """

        # FIXME: make it model independent

        if gf is None:

            gf = self.gf_low

        assert gf > 0 and gf <= 1.5



        data = zip(data.tissues, self.N2_A, self.N2_B, self.HE_A, self.HE_B)

        return tuple(

            eq_gf_limit(gf, p_n2, p_he, n2_a, n2_b, he_a, he_b)

            for (p_n2, p_he), n2_a, n2_b, he_a, he_b in data

        )





class ZH_L16B_GF(ZH_L16_GF): # source: gfdeco.f by Baker

    """

    ZH-L16B-GF decompression model.

    """

    N2_A = (

        1.1696, 1.0000, 0.8618, 0.7562, 0.6667, 0.5600, 0.4947, 0.4500,

        0.4187, 0.3798, 0.3497, 0.3223, 0.2850, 0.2737, 0.2523, 0.2327,

    )

    N2_B = (

        0.5578, 0.6514, 0.7222, 0.7825, 0.8126, 0.8434, 0.8693, 0.8910,

        0.9092, 0.9222, 0.9319, 0.9403, 0.9477, 0.9544, 0.9602, 0.9653,

    )

    HE_A = (

        1.6189, 1.3830, 1.1919, 1.0458, 0.9220, 0.8205, 0.7305, 0.6502,

        0.5950, 0.5545, 0.5333, 0.5189, 0.5181, 0.5176, 0.5172, 0.5119,

    )

    HE_B = (

        0.4770, 0.5747, 0.6527, 0.7223, 0.7582, 0.7957, 0.8279, 0.8553,

        0.8757, 0.8903, 0.8997, 0.9073, 0.9122, 0.9171, 0.9217, 0.9267,

    )

    N2_HALF_LIFE = (

        5.0, 8.0, 12.5, 18.5, 27.0, 38.3, 54.3, 77.0, 109.0,

        146.0, 187.0, 239.0, 305.0, 390.0, 498.0, 635.0,

    )

    HE_HALF_LIFE = (

        1.88, 3.02, 4.72, 6.99, 10.21, 14.48, 20.53, 29.11,

        41.20, 55.19, 70.69, 90.34, 115.29, 147.42, 188.24, 240.03,

    )







class ZH_L16C_GF(ZH_L16_GF): # source: ostc firmware code

    """

    ZH-L16C-GF decompression model.

    """

    N2_A = (

        1.2599, 1.0000, 0.8618, 0.7562, 0.6200, 0.5043, 0.4410, 0.4000,

        0.3750, 0.3500, 0.3295, 0.3065, 0.2835, 0.2610, 0.2480, 0.2327,

    )

    N2_B = (

        0.5050, 0.6514, 0.7222, 0.7825, 0.8126, 0.8434, 0.8693, 0.8910,

        0.9092, 0.9222, 0.9319, 0.9403, 0.9477, 0.9544, 0.9602, 0.9653,

    )

    HE_A = (

        1.7424, 1.3830, 1.1919, 1.0458, 0.9220, 0.8205, 0.7305, 0.6502,

        0.5950, 0.5545, 0.5333, 0.5189, 0.5181, 0.5176, 0.5172, 0.5119,

    )

    HE_B = (

        0.4245, 0.5747, 0.6527, 0.7223, 0.7582, 0.7957, 0.8279, 0.8553,

        0.8757, 0.8903, 0.8997, 0.9073, 0.9122, 0.9171, 0.9217, 0.9267,

    )

    N2_HALF_LIFE = (

        4.0, 8.0, 12.5, 18.5, 27.0, 38.3, 54.3, 77.0, 109.0,

        146.0, 187.0, 239.0, 305.0, 390.0, 498.0, 635.0,

    )

    HE_HALF_LIFE = (

        1.51, 3.02, 4.72, 6.99, 10.21, 14.48, 20.53, 29.11, 41.20,

        55.19, 70.69, 90.34, 115.29, 147.42, 188.24, 240.03,

    )



