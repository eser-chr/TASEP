import numpy as np


class Analytical:
    @staticmethod
    def calc_P(kon: float, koff: float, kstep: float, grid: np.ndarray) -> np.ndarray:
        """
        Calculate probabilities using formula eq.10

        Parameters:
            kon (float): On-rate.
            koff (float): Off-rate.
            kstep (float): Step rate.
            grid (numpy.ndarray): The input grid.

        Returns:
            numpy.ndarray: Probabilities.
        """
        L = len(grid)
        n = np.arange(0, L)
        return (kon / koff) * (1 - np.exp(-koff * n / kstep))

    @staticmethod
    def num_kinesins(L: int, kon: float, koff: float, kstep: float) -> np.float64:
        """
        Calculate the number of kinesins, as a discrete sum of eq.10 over the sides.

        Parameters:
            L (int): Length of the grid.
            kon (float): On-rate.
            koff (float): Off-rate.
            kstep (float): Step rate.

        Returns:
            numpy.float64: Number of kinesins.
        """
        n = np.arange(0, L)
        return np.sum((kon / koff) * (1 - np.exp(-koff * n / kstep)))