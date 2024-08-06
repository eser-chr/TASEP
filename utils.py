import matplotlib.pyplot as plt
import numpy as np
from joblib import Parallel, delayed
import sys
import cooperative_tasep_lib as tasep



class Plot:
    @staticmethod        
    def nice_figure(data):
        fig, ax = plt.subplots(1,1, figsize=(10,10))
        ax.imshow(data)
        ax.set(title = "density-gram", xlabel = "microtrubule (sides)", ylabel="time (A.U)")

    @staticmethod
    def total_density_time(data):
        S = np.mean(data, axis=1)
        fig, ax = plt.subplots(1,1, figsize=(10,10))
        ax.plot(S)
        ax.set(title = "Kins vs Time", xlabel = "time (A.U)", ylabel="density")


class CollectionPlot:
    @staticmethod
    def total_density_time(data):
        S = np.mean(data, axis=1)
        S = np.mean(S, axis=1)
        fig, ax = plt.subplots(1,1, figsize=(10,10))
        ax.plot(S)
        ax.set(title = "Kins vs Time", xlabel = "time (A.U)", ylabel="density")



class MultipleExecution:
    @staticmethod
    def same_conf_parallel(L, T, kon, koff, kstep, kq, q, num=100):
        results = Parallel(n_jobs=-1)(
            delayed(tasep.ssim)(L, T, kon, koff, kstep, kq, q)
            for _ in range(num)
        )
        
        return [res[0] for res in results]
    
    def same_conf_serial(L, T, kon, koff, kstep, kq, q, num=100):
        RES = np.zeros((10*T+1, L , num))

        for _ in range(num):
            data, times = tasep.sim(L, T, kon, koff, kstep, kq, q)
            sys.stdout.flush()
            RES[:, :, _] = data

        return RES

    
    @staticmethod
    def sweep_over():
        pass