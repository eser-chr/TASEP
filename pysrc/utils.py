import matplotlib.pyplot as plt
import numpy as np


class DataMismatchError(Exception):
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)







class CollectionAnalysis:
    @staticmethod
    def total_density_time(data):
        S = np.mean(data, axis=1)
        S = np.mean(S, axis=1)
        return S
    
    @staticmethod
    def get_val_at_time(data, time):
        return CollectionAnalysis.total_density_time(data)[time]
    
    @staticmethod
    def normalize_at_time(data, time):
        timeseries = CollectionAnalysis.total_density_time(data)
        return timeseries/timeseries[time]


class Theoretical:
    class evolution:
        def non_coop(dt, tmax, kon, koff):
            f = kon/(kon+koff)
            x = np.arange(0, tmax, dt)
            y = f*(1-np.exp(-(kon+koff)*x))
            return x, y
        
        def coop(dt, tmax, kon, koff, q, a):
            x = np.arange(0, tmax, dt)
            y = np.zeros_like(x)
            Q = q-1
            for i in range(len(x)-1):
                y[i+1] = y[i] + dt * (kon*(1-y[i]) + Q*kon*a[i] - koff*y[i])

            return x, y




class CollectionPlots:
    @staticmethod
    def total_density_time(data, labels=[], size = (10,10)):
        if len(data)!=len(labels):
            raise DataMismatchError("length of data is not the same as len of labels.")
        
        fig, ax = plt.subplots(1,1, figsize=size)
        for i, line in enumerate(data):
            S = np.mean(line, axis=1)
            S = np.mean(S, axis = 1)
            ax.plot(S, label= labels[i])

        ax.set(title = f"Kins vs Time", xlabel = "time (A.U)", ylabel="density")
        ax.legend()
        ax.grid()
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
    
            



