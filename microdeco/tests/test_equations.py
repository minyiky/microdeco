import unittest

from microdeco.equations import schreiner, buhlmann

class TestSchreinerFunction(unittest.TestCase):

    def test_descent(self):
        self.assertAlmostEqual(schreiner(p_alv=0.637364, p_t=0.74065446, k=0.138629, t=1.5, R=1.36), 0.919397, 5)

    def test_constant_depth(self):
        self.assertAlmostEqual(schreiner(p_alv=2.677364, p_t=0.919397, k=0.138629, t=20, R=0), 2.567490, 5)

    def test_ascent(self):
        self.assertAlmostEqual(schreiner(p_alv=2.677364, p_t=2.567490, k=0.138629, t=2, R=-0.68), 2.421840, 5)

class TestBuhlmannFunction(unittest.TestCase):

    def test_surface(self):
        self.assertAlmostEqual(buhlmann(0.3, 0.74065446, 0, 1.1696, 0, 0.5578, 0), 0.314886, 5)

    def test_descent(self):
        self.assertAlmostEqual(buhlmann(0.3, 0.919397, 0, 1.1696, 0, 0.5578, 0), 0.4592862, 5)

    def test_constant_depth(self):
        self.assertAlmostEqual(buhlmann(0.3, 2.567490, 0, 1.1696, 0, 0.5578, 0), 1.7907266, 5)

    def test_ascent(self):
        self.assertAlmostEqual(buhlmann(0.3, 2.421840, 0, 1.1696, 0, 0.5578, 0), 1.6730607, 5)


if __name__ == '__main__':
    unittest.main()
