import pickle
import numpy as np
from typing import Callable

def write_sim(path:str, data:np.ndarray, times:np.ndarray)->bool:
    """
    Writes all the data in pickle form

    Args:
        path: The path of the file 
        data: The data of the simulation
        times: the time points of the simulation
    
    Returns:
        True: if the writing was succesful

    """

    with open(path, 'wb') as f:
        pickle.dump({"data":data, "times":times}, f)

    return True


def moving_func(grid: np.ndarray, window_size: int, 
                func: Callable[[np.ndarray], np.ndarray]) -> np.ndarray:
    """
    Calculate a moving function over a grid.

    Parameters:
        grid (numpy.ndarray): The input grid.
        window_size (int): The size of the moving window.
        func (callable): The function to be applied at each slice. Ideally its a numpy function such as np.mean

    Returns:
        numpy.ndarray: The result of applying the function over the moving window.
    """
    moving = np.zeros_like(grid)
    for i in range(len(grid)):
        start_index = max(0, i - window_size // 2)
        end_index = min(len(grid), i + window_size // 2 + 1)
        moving[i] = func(grid[start_index:end_index])
    return moving