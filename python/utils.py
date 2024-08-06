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

def secure_ssim(L, T, kon, koff, kstep, kq, q):
    data, activation, nn, times, res, dts = tasep.ssim(L, T, kon, koff, kstep, kq, q)
    sys.stdout.flush()
    return data, activation


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
    
            



class MultipleExecutionSimple:    
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

class MultipleExecutionDetailed:    
    @staticmethod
    def same_conf_parallel(L, T, kon, koff, kstep, kq, q, num=100):
        results = Parallel(n_jobs=10, prefer="threads")(
            delayed(secure_ssim)(L, T, kon, koff, kstep, kq, q)
            for _ in range(num)
        )
        data = [res[0] for res in results]
        activation = [res[1] for res in results]
        return (np.stack(data, axis=-1), np.stack(activation, axis=-1))
        # return np.stack(results, axis=-1)

    
    # def same_conf_serial(L, T, kon, koff, kstep, kq, q, num=100):
    #     RES = np.zeros((10*T+1, L , num))
    #     for _ in range(num):
    #         data = secure_ssim(L, T, kon, koff, kstep, kq, q)
    #         sys.stdout.flush()
    #         RES[:, :, _] = data

    #     return RES

    
    # @staticmethod
    # def sweep_over():
    #     pass