from binding import _line_sim, _cyclic_sim
import numpy as np
from typing import Tuple

def line_sim(L: int, T: float, kon: float, koff: float, kstep: float, kq: float, q: float) -> Tuple[np.ndarray,np.ndarray,np.ndarray,np.ndarray ]:
    """
    Simulate a line system.

    Parameters:
    - L (int): Length of the line.
    - T (float): Total simulation time.
    - kon (float): Rate constant for kinesin binding.
    - koff (float): Rate constant for kinesin unbinding.
    - kstep (float): Rate constant for kinesin stepping.
    - kq (float): Rate constant for external force.
    - q (float): Magnitude of external force.

    Returns:
    - tuple: A tuple containing four NumPy arrays:
        - data (np.ndarray): A 2D array representing the state of the line system over time.
        - times (np.ndarray): An array containing the times at which the system state was recorded.
        - res (np.ndarray): An array containing the results of the simulation.
        - dt (np.ndarray): An array containing the time steps used in the simulation.
    """
    # Simulate the line system
    data, times, res, dt = _line_sim(L, T, kon, koff, kstep, kq, q)

    # Convert lists to NumPy arrays
    data = np.array(data)
    times = np.array(times)
    res = np.array(res)
    dt = np.array(dt)

    return data, times, res, dt

def cyclic_sim(L: int, T: float, kon: float, koff: float, kstep: float, kq: float, q: float, initial_density:int=0) -> Tuple[np.ndarray,np.ndarray,np.ndarray,np.ndarray ]:
    """
    Simulate a system (S1 topology/periodic bdrs).

    Parameters:
    - L (int): Length of the line.
    - T (float): Total simulation time.
    - kon (float): Rate constant for kinesin binding.
    - koff (float): Rate constant for kinesin unbinding.
    - kstep (float): Rate constant for kinesin stepping.
    - kq (float): Rate constant for external force.
    - q (float): Magnitude of external force.
    - initial_density (int): The percentage of the initial side occupation

    Returns:
    - tuple: A tuple containing four NumPy arrays:
        - data (np.ndarray): A 2D array representing the state of the line system over time.
        - times (np.ndarray): An array containing the times at which the system state was recorded.
        - res (np.ndarray): An array containing the results of the simulation.
        - dt (np.ndarray): An array containing the time steps used in the simulation.
    """
    # Simulate the line system
    data, times, res, dt = _cyclic_sim(L, T, kon, koff, kstep, kq, q, initial_density)

    # Convert lists to NumPy arrays
    data = np.array(data)
    times = np.array(times)
    res = np.array(res)
    dt = np.array(dt)

    return data, times, res, dt
