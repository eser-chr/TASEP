import numpy as np
import matplotlib.pyplot as plt


"""
Plot Class

The Plot class provides static methods to generate specific types of plots using matplotlib.

Methods
-------
nice_figure(data)
    Generates a heatmap of the given data.
    
    Parameters:
    - data (numpy array): The data to be plotted as a heatmap.
    
    Returns:
    - matplotlib.figure.Figure: The figure object containing the heatmap.
    
    Example:
    --------
    import numpy as np
    import matplotlib.pyplot as plt
    
    data = np.random.random((100, 100))
    fig = Plot.nice_figure(data)
    plt.show()
    
    Notes:
    ------
    - The heatmap represents the density of data.
    - The x-axis is labeled "microtubule (sides)".
    - The y-axis is labeled "time (A.U)".
    - The title of the plot is "density-gram".

total_density_time(data)
    Generates a line plot showing the average density over time.
    
    Parameters:
    - data (numpy array): The data to be plotted. The mean is calculated along the second axis (axis=1).
    
    Returns:
    - matplotlib.figure.Figure: The figure object containing the line plot.
    
    Example:
    --------
    import numpy as np
    import matplotlib.pyplot as plt
    
    data = np.random.random((100, 100))
    fig = Plot.total_density_time(data)
    plt.show()
    
    Notes:
    ------
    - The x-axis is labeled "time (A.U)".
    - The y-axis is labeled "density".
    - The title of the plot is "Kins vs Time".
"""

class Plot:
    @staticmethod        
    def nice_figure(data):
        fig, ax = plt.subplots(1,1, figsize=(10,10))
        ax.imshow(data)
        ax.set(title="density-gram", xlabel="microtrubule (sides)", ylabel="time (A.U)")
        return fig

    @staticmethod
    def total_density_time(data):
        S = np.mean(data, axis=1)
        fig, ax = plt.subplots(1,1, figsize=(10,10))
        ax.plot(S)
        ax.set(title="Kins vs Time", xlabel="time (A.U)", ylabel="density")
        return fig

"""
CollectionPlot Class

The CollectionPlot class provides a static method to generate a line plot showing the average density over time 
across multiple simulations.

Methods
-------
total_density_time(data, size=(10,10))
    Generates a line plot showing the average density over time across multiple simulations.
    
    Parameters:
    - data (numpy array): The data to be plotted. The mean is calculated along the second axis (axis=1) 
      and the third axis (axis=2).
    - size (tuple): The size of the figure (default is (10,10)).
    
    Returns:
    - matplotlib.figure.Figure: The figure object containing the line plot.
    
    Example:
    --------
    import numpy as np
    import matplotlib.pyplot as plt
    
    data = np.random.random((100, 100, 10))
    fig = CollectionPlot.total_density_time(data)
    plt.show()
    
    Notes:
    ------
    - The x-axis is labeled "time (A.U)".
    - The y-axis is labeled "density".
    - The title of the plot includes the average number of simulations.
"""

class CollectionPlot:
    @staticmethod
    def total_density_time(data, size=(10,10)):
        num_of_sims = data.shape[2]
        S = np.mean(data, axis=1)
        S = np.mean(S, axis=1)
        fig, ax = plt.subplots(1,1, figsize=size)
        ax.plot(S)
        ax.set(title=f"Kins vs Time (avg:{num_of_sims} sims)", xlabel="time (A.U)", ylabel="density")
        return fig
    
    @staticmethod
    def timeseries(data, labels, size=(10,10)):
        if len(data)!=len(labels):
            raise DataMismatchError("length of data is not the same as len of labels.")

        fig, ax = plt.subplots(1,1, figsize=size)
        for i, line in enumerate(data):
            ax.plot(line, label= labels[i])

        ax.set(title = f"Kins vs Time", xlabel = "time (A.U)", ylabel="density")
        ax.legend()
        ax.grid()
        return fig








