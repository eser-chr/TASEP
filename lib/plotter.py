import numpy as np
import matplotlib.pyplot as plt
from lib.utils import moving_func
from lib.analytical import Analytical

def frame(grid: np.ndarray, kon: float, koff: float, kstep: float, sims: int, window: int = 1) -> None:
    """
    Plot simulation and analytical data.

    Parameters:
        grid (numpy.ndarray): The input grid.
        kon (float): On-rate.
        koff (float): Off-rate.
        kstep (float): Step rate.
        sims (int): Number of simulations.
        window (int, optional): Size of the moving window. Defaults to 1.
    """
    avg = moving_func(grid, window, np.mean)
    std = moving_func(grid, window, np.std)

    theo3 = Analytical.calc_P(kon, koff, kstep, grid)   

    fig, ax = plt.subplots(1, 1, figsize=(10, 5))
    ax.plot(avg, label='simulation')
    ax.plot(theo3, c='r', linestyle='-.', label='analytical')

    ax.fill_between(range(len(avg)), avg - std, avg + std, color='lightblue', alpha=0.5)
    ax.set(title=f'Total num of kinesins={int(np.sum(grid))} \n Average over {sims} sims', ylim=[0, 4])
    ax.grid()
    ax.legend(loc=(0.805, 0.125))
    plt.show()