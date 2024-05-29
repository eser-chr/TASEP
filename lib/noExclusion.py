import numpy as np
from numba import njit, prange
from typing import Tuple

@njit
def calc_propensities(grid: np.ndarray[np.any, np.dtype[np.float64]], 
                      kon: float, 
                      koff: float, 
                      kstep: float) -> np.ndarray[np.float64]:
    """
    Calculate propensities for different events.

    Args:
        grid (np.ndarray): 1D array representing the grid state.
        kon (float): Association rate constant.
        koff (float): Dissociation rate constant.
        kstep (float): Step rate constant.

    Returns:
        np.ndarray: Array containing propensities for association, dissociation, and step events.
    """
    n_on = len(grid)
    n_off = np.sum(grid)
    n_step = np.sum(grid)

    aon = n_on * kon
    aoff = n_off * koff
    astep = n_step * kstep

    return np.array([aon, aoff, astep])

@njit
def bind_kinesin(grid: np.ndarray[np.any, np.dtype[np.float64]]) -> None:
    """
    Bind a kinesin to a random site on the grid.

    Args:
        grid (np.ndarray): 1D array representing the grid state.
    """
    side = np.random.randint(len(grid))
    grid[side] += 1

@njit
def unbind_kinesin(grid: np.ndarray[np.any, np.dtype[np.float64]]) -> None:
    """
    Unbind a kinesin from a random site on the grid.

    Args:
        grid (np.ndarray): 1D array representing the grid state.
    """
    S = np.cumsum(grid)
    r = np.random.randint(low=1, high=int(S[-1]) + 1)
    side = np.argmax(S >= r)

    grid[side] -= 1

@njit
def move_kinesin_w_fall(grid: np.ndarray[np.any, np.dtype[np.float64]]) -> None:
    """
    Move a kinesin along the grid with a possibility of falling.

    Args:
        grid (np.ndarray): 1D array representing the grid state.
    """
    S = np.cumsum(grid)
    r = np.random.randint(low=1, high=int(S[-1]) + 1)
    side = np.argmax(S >= r)
    
    grid[side] -= 1
    if side != len(grid) - 1:
        grid[side + 1] += 1

@njit
def step(grid: np.ndarray[np.any, np.dtype[np.float64]], 
         kon: float, 
         koff: float, 
         kstep: float) -> Tuple[float, int]:
    """
    Perform a step in the simulation.

    Args:
        grid (np.ndarray): 1D array representing the grid state.
        kon (float): Association rate constant.
        koff (float): Dissociation rate constant.
        kstep (float): Step rate constant.

    Returns:
        Tuple[float, int]: A tuple containing the time step and the index of the event.
    """
    r1 = np.random.uniform()
    r2 = np.random.uniform()

    A = calc_propensities(grid, kon, koff, kstep)
    R_tot = np.sum(A)
    
    A_normalised = np.cumsum(A) / R_tot

    dt = (1 / R_tot) * np.log(1 / r1)
    idx = np.argwhere(A_normalised > r2)[0][0]

    return dt, int(idx)


@njit
def simulation(L: int, T: float, kon: float, koff: float, kstep: float) -> Tuple[np.ndarray, np.ndarray]:
    """
    Run a simulation of kinesin movement on a grid.

    Args:
        L (int): Length of the grid.
        T (float): Total simulation time.
        kon (float): Association rate constant.
        koff (float): Dissociation rate constant.
        kstep (float): Step rate constant.

    Returns:
        Tuple[np.ndarray, np.ndarray]: A tuple containing data array and time array.
    """
    time = 0.0
    block = 0
    period = 0.1
    next_write_time = period
    blocks = int(T / period)

    grid = np.zeros(L, dtype=np.float64)

    DATA = np.zeros((blocks, L), dtype=np.float64)
    TIMES = np.zeros(blocks, dtype=np.float64)

    while time < T:
        dt, idx = step(grid, kon, koff, kstep)
        time += dt      
        if idx == 0:
            bind_kinesin(grid)
        elif idx == 1:
            unbind_kinesin(grid)
        else:
            move_kinesin_w_fall(grid)
            
        if next_write_time < time:
            DATA[block, :] = grid
            TIMES[block] = time
            block += 1
            next_write_time += period

    return DATA[:block], TIMES[:block]


@njit(parallel=True)
def mul_sim_last(num:int, L:int, T:int, kon:float, koff:float, 
                 kstep:float)->np.ndarray:
    """ 
    Runs multiple of simulations and returns an array of only the last frame

    Args: 
        num: Number of simulations
        L (int): Length of the grid.
        T (float): Total simulation time.
        kon (float): Association rate constant.
        koff (float): Dissociation rate constant.
        kstep (float): Step rate constant.

    Returns:
        Numpy array [Number of sims, Length of grid]
    """

    LAST_FRAMES = np.zeros((num, L))
    for i in prange(num):
        DATA, TIMES = simulation(L, T, kon, koff, kstep)
        LAST_FRAMES[i,:]=DATA[-1,:]
    return LAST_FRAMES


@njit(parallel=True)
def mul_sim_all(num:int, L:int, T:int, kon:float, koff:float, kstep:float
                  )->tuple[np.ndarray, np.ndarray]:
    """ 
    Runs multiple of simulations and returns an array of only the last frame

    Args: 
        num: Number of simulations
        L (int): Length of the grid.
        T (float): Total simulation time.
        kon (float): Association rate constant.
        koff (float): Dissociation rate constant.
        kstep (float): Step rate constant.

    Returns:
        Numpy array [Number of sims, Length of grid]
    """

    DATA = np.zeros((num, int(10*T), L))
    TIMES = np.zeros((num, int(10*T))) # Because a period is 0.1. So we have at most 10*T
    for i in prange(num):
        data, times = simulation(L, T, kon, koff, kstep)
        max_T, max_L = data.shape
        DATA[i,:max_T, :]=data[:max_T,:]
        TIMES[i, :max_T] = times
    return DATA, TIMES