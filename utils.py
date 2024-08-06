import matplotlib.pyplot as plt
import numpy as np
from joblib import Parallel, delayed
import sys
import cooperative_tasep_lib as tasep

class DataMismatchError(Exception):
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)


def secure_sim(L, T, kon, koff, kstep, kq, q):
    data, times = tasep.sim(L, T, kon, koff, kstep, kq, q)
    sys.stdout.flush()
    return data


class Plot:
    @staticmethod        
    def nice_figure(data):
        fig, ax = plt.subplots(1,1, figsize=(10,10))
        ax.imshow(data)
        ax.set(title = "density-gram", xlabel = "microtrubule (sides)", ylabel="time (A.U)")
        return fig

    @staticmethod
    def total_density_time(data):
        S = np.mean(data, axis=1)
        fig, ax = plt.subplots(1,1, figsize=(10,10))
        ax.plot(S)
        ax.set(title = "Kins vs Time", xlabel = "time (A.U)", ylabel="density")
        return fig


class CollectionPlot:
    @staticmethod
    def total_density_time(data, size = (10,10)):
        num_of_sims = data.shape[2]
        S = np.mean(data, axis=1)
        S = np.mean(S, axis=1)
        fig, ax = plt.subplots(1,1, figsize=size)
        ax.plot(S)
        ax.set(title = f"Kins vs Time (avg:{num_of_sims} sims)", xlabel = "time (A.U)", ylabel="density")
        return fig
    


    

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
        return fig
    
            



class MultipleExecution:    
    @staticmethod
    def same_conf_parallel(L, T, kon, koff, kstep, kq, q, num=100):
        results = Parallel(n_jobs=10, prefer="threads")(
            delayed(secure_sim)(L, T, kon, koff, kstep, kq, q)
            for _ in range(num)
        )
        
        return np.stack(results, axis=-1)
    
    def same_conf_serial(L, T, kon, koff, kstep, kq, q, num=100):
        RES = np.zeros((10*T+1, L , num))
        for _ in range(num):
            data = secure_sim(L, T, kon, koff, kstep, kq, q)
            sys.stdout.flush()
            RES[:, :, _] = data

        return RES

    
    @staticmethod
    def sweep_over():
        pass