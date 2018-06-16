import unittest
import numpy.testing
import asymk

import numpy as np

class TestAsymK(unittest.TestCase):
    """
    Test class of asymk.py.
    """

    def test_DeriveGraph(self):
        """test method for tashizan
        """
        print('Test DeriveGraph')
        D = np.matrix([[0, 3], [2, 0]])
        R = 2
        expected = np.matrix([[0, 0], [1, 0]])
        actual = asymk.DeriveGraph(D, R)
        print('expected=', expected)
        print('actual  =', actual)
        numpy.testing.assert_equal(expected, actual)

    def test_Gamma(self):
        print('Test Gamma')
        Gr = np.matrix([[0, 1, 0], [0, 0, 1], [0, 0, 0]], dtype=int)

        # Gamma Plus
        P = np.array([1, 0, 0])
        R1 = 1
        expectedA = np.matrix([[1, 1, 0]])
        actualA = asymk.Gamma(Gr, P, R1, True)
        numpy.testing.assert_equal(expectedA, actualA)

        R2 = 2
        expectedB = np.matrix([[1, 1, 1]])
        actualB = asymk.Gamma(Gr, P, R2, True)
        numpy.testing.assert_equal(expectedB, actualB)

        # Gamma Minus
        Pc = np.array([0, 0, 1])
        expectedC = np.matrix([[0, 1, 1]])
        actualC = asymk.Gamma(Gr, Pc, R1, False)
        numpy.testing.assert_equal(expectedC, actualC)

        expectedD = np.matrix([[1, 1, 1]])
        actualD = asymk.Gamma(Gr, Pc, R2, False)
        numpy.testing.assert_equal(expectedD, actualD)


    def test_CCV(self):
        print('Test CCV')
        Gr = np.matrix([[0, 1, 0], [1, 0, 1], [0, 0, 0]], dtype=int)
        A  = np.array([1, 0, 0])
        expectedA = np.array([1, 0, 0])
        actualA = asymk.CCV(Gr, A)
        numpy.testing.assert_equal(expectedA, actualA)

if __name__ == "__main__":
    unittest.main()
