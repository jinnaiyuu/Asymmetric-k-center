import unittest
import numpy.testing
import asymk

import numpy as np

class TestAsymK(unittest.TestCase):
    """test class of asymk.py
    """

    def test_DeriveGraph(self):
        """test method for tashizan
        """
        print('Test DeriveGraph')
        D = np.matrix([[0, 3], [2, 0]])
        R = 2
        expected = np.matrix([[1, 0], [1, 1]])
        actual = asymk.DeriveGraph(D, R)
        print('expected=', expected)
        print('actual  =', actual)
        numpy.testing.assert_equal(expected, actual)

    def test_GammaPlus(self):
        print('Test GammaPlus')
        Gr = np.matrix([[0, 1, 0], [0, 0, 1], [0, 0, 0]], dtype=int)
        P = np.array([1, 0, 0])
        R1 = 1
        expectedA = np.matrix([[1, 1, 0]])
        actualA = asymk.GammaPlus(Gr, P, R1)
        numpy.testing.assert_equal(expectedA, actualA)

        R2 = 2
        expectedB = np.matrix([[1, 1, 1]])
        actualB = asymk.GammaPlus(Gr, P, R2)
        numpy.testing.assert_equal(expectedB, actualB)

    def test_CCV(self):
        print('Test CCV')

    def test_Reduce(Gr):
        print('Test Reduce')

    def test_Cover(G, C):
        print('Test Cover')

    def test_LP(G, A):
        print('Test LP')
        
if __name__ == "__main__":
    unittest.main()
