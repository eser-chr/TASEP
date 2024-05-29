from numba import njit, prange
import numpy as np
from typing import Tuple

@njit
def activation_to_factors(L: int, q: float, grid_activated: np.ndarray) -> np.ndarray:
    """
    Calculate factors based on activation.

    Args:
        L (int): Length of the grid.
        q (float): Activation factor.
        grid_activated (np.ndarray): 1D array representing the activated grid.

    Returns:
        np.ndarray: Array of factors.
    """
    return np.ones(L) + (q - 1) * grid_activated

@njit
def calc_propensities(L: int, grid: np.ndarray, grid_activated: np.ndarray, kon: float, koff: float, kstep: float, kq: float, q: float) -> np.ndarray:
    """
    Calculate propensities for different events.

    Args:
        L (int): Length of the grid.
        grid (np.ndarray): 1D array representing the grid state.
        grid_activated (np.ndarray): 1D array representing the activated grid state.
        kon (float): Association rate constant.
        koff (float): Dissociation rate constant.
        kstep (float): Step rate constant.
        kq (float): Activation rate constant.
        q (float): Activation factor.

    Returns:
        np.ndarray: Array containing propensities for association, dissociation, step, and activation events.
    """
    n_on = np.sum(activation_to_factors(L, q, grid_activated))
    n_off = np.sum(grid)
    n_step = np.sum(grid)
    n_q = len(np.where((grid_activated == 1) & (grid == 0))[0])

    aon = n_on * kon
    aoff = n_off * koff
    astep = n_step * kstep
    aq = n_q * kq

    return np.array([aon, aoff, astep, aq])

@njit
def bind_kinesin(L: int, q: float, grid: np.ndarray, grid_activated: np.ndarray) -> None:
    """
    Bind a kinesin to a random site on the grid.

    Args:
        L (int): Length of the grid.
        q (float): Activation factor.
        grid (np.ndarray): 1D array representing the grid state.
        grid_activated (np.ndarray): 1D array representing the activated grid state.
    """
    S = np.cumsum(activation_to_factors(L, q, grid_activated))
    r = S[-1] * np.random.random()
    side = np.argmax(S >= r)
    grid[side] += 1

@njit
def unbind_kinesin(grid: np.ndarray) -> None:
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
def move_kinesin_w_fall(grid: np.ndarray) -> None:
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
def deactivate_tubule(grid: np.ndarray, grid_activated: np.ndarray) -> None:
    """
    Deactivate a tubule on the grid.

    Args:
        grid (np.ndarray): 1D array representing the grid state.
        grid_activated (np.ndarray): 1D array representing the activated grid state.
    """
    sides = np.where((grid_activated == 1) & (grid == 0))[0]
    S = len(sides)
    r = np.random.randint(S)
    side = sides[r]
    grid_activated[side] = 0

@njit
def activate_tubule(grid: np.ndarray, grid_activated: np.ndarray) -> np.ndarray:
    """
    Activate a tubule on the grid.

    Args:
        grid (np.ndarray): 1D array representing the grid state.
        grid_activated (np.ndarray): 1D array representing the activated grid state.

    Returns:
        np.ndarray: 1D array representing the updated activated grid state.
    """
    sides = np.where(grid != 0)[0]
    grid_activated[sides] = 1
    return grid_activated

@njit
def step(L: int, grid: np.ndarray, grid_activated: np.ndarray, kon: float, koff: float, kstep: float, kq: float, q: float) -> Tuple[float, int]:
    """
    Perform a step in the simulation.

    Args:
        L (int): Length of the grid.
        grid (np.ndarray): 1D array representing the grid state.
        grid_activated (np.ndarray): 1D array representing the activated grid state.
        kon (float): Association rate constant.
        koff (float): Dissociation rate constant.
        kstep (float): Step rate constant.
        kq (float): Activation rate constant.
        q (float): Activation factor.

    Returns:
        Tuple[float, int]: A tuple containing the time step and the index of the event.
    """
    r1 = np.random.uniform()
    r2 = np.random.uniform()

    A = calc_propensities(L, grid, grid_activated, kon, koff, kstep, kq, q)
    R_tot = np.sum(A)
    
    A_normalised = np.cumsum(A) / R_tot

    dt = (1 / R_tot) * np.log(1 / r1)
    idx = np.argwhere(A_normalised > r2)[0][0]
    return dt, int(idx)



@njit
def simulation(L: int, T: float, kon: float, koff: float, kstep: float, kq: float, q: float) -> Tuple[np.ndarray, np.ndarray]:
    """
    Run a simulation of kinesin movement on a grid.

    Args:
        L (int): Length of the grid.
        T (float): Total simulation time.
        kon (float): Association rate constant.
        koff (float): Dissociation rate constant.
        kstep (float): Step rate constant.
        kq (float): Activation rate constant.
        q (float): Activation factor.

    Returns:
        Tuple[np.ndarray, np.ndarray]: A tuple containing data array and time array.
    """
    time = 0.0
    block = 0
    period = 0.1
    next_write_time = period
    blocks = int(T / period)

    grid = np.ones(L, dtype=float)
    grid_activated = np.zeros(L, dtype=float)

    DATA = np.zeros((blocks, L))
    TIMES = np.zeros(blocks)

    while time < T:
        dt, idx = step(L, grid, grid_activated, kon, koff, kstep, kq, q)
        time += dt      

        # Action
        if idx == 0:
            bind_kinesin(L, q, grid, grid_activated)
        elif idx == 1:
            unbind_kinesin(grid)
        elif idx == 2:
            move_kinesin_w_fall(grid)
        else:
            deactivate_tubule(grid, grid_activated)

        grid_activated = activate_tubule(grid, grid_activated)
        
        # Write    
        if next_write_time < time:
            DATA[block, :] = grid
            TIMES[block] = time
            block += 1
            next_write_time += period

    return DATA[:block], TIMES[:block]


@njit(parallel=True)
def mul_sim_last(num:int, L:int, T:int, kon:float, koff:float, kstep:float,
                  kq: float, q: float)->np.ndarray:
    """ 
    Runs multiple of simulations and returns an array of only the last frame

    Args: 
        num: Number of simulations
        L (int): Length of the grid.
        T (float): Total simulation time.
        kon (float): Association rate constant.
        koff (float): Dissociation rate constant.
        kstep (float): Step rate constant.
        kq (float): Activation rate constant.
        q (float): Activation factor.

    Returns:
        Numpy array [Number of sims, Length of grid]
    """

    LAST_FRAMES = np.zeros((num, L))
    for i in prange(num):
        DATA, TIMES = simulation(L, T, kon, koff, kstep, kq, q)
        LAST_FRAMES[i,:]=DATA[-1,:]
    return LAST_FRAMES


@njit(parallel=True)
def mul_sim_all(num:int, L:int, T:int, kon:float, koff:float, kstep:float,
                  kq: float, q: float)->tuple[np.ndarray, np.ndarray]:
    """ 
    Runs multiple of simulations and returns an array of only the last frame

    Args: 
        num: Number of simulations
        L (int): Length of the grid.
        T (float): Total simulation time.
        kon (float): Association rate constant.
        koff (float): Dissociation rate constant.
        kstep (float): Step rate constant.
        kq (float): Activation rate constant.
        q (float): Activation factor.

    Returns:
        Numpy array [Number of sims, Length of grid]
    """

    DATA = np.zeros((num, int(10*T), L))
    TIMES = np.zeros((num, int(10*T))) # Because a period is 0.1. So we have at most 10*T
    for i in prange(num):
        data, times = simulation(L, T, kon, koff, kstep, kq, q)
        max_T, max_L = data.shape
        DATA[i,:max_T, :]=data[:max_T,:]
        TIMES[i, :max_T] = times
    return DATA, TIMES



