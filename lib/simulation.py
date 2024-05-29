import os
import pickle
import pandas as pd
from lib.parameters import Parameters
from lib.paths import Paths

class Simulation:
    """
    A class to handle simulations and related operations.

    Attributes:
        parameters (Parameters): The parameters object for the simulation.
        write (bool): Flag indicating whether to write the simulation results.
        exclusion (bool): Flag indicating whether exclusion is enabled.
        _sim (function): Function reference for simulation.
        _mul_sim (function): Function reference for multiple simulations with last frame output.
        _mul_sim_all (function): Function reference for multiple simulations with all frames output.
        options (tuple): Simulation options as a tuple.
    """

    def __init__(self, parameters: Parameters, exclusion: bool, write: bool):
        """
        Initializes the Simulation object.

        Args:
            parameters (Parameters): The parameters object for the simulation.
            exclusion (bool): Flag indicating whether exclusion is enabled.
            write (bool): Flag indicating whether to write the simulation results.
        """
        if parameters.no_debug_check():
            self.parameters = parameters
        else: 
            parameters.debug_check()

        self.write = write
        self.exclusion = exclusion
        self._sim = None
        self._mul_sim = None
        self._mul_sim_all = None

    def configure(self):  
        """
        Configures the simulation based on parameters and exclusion flag.
        """
        if self.parameters.q == 1 and self.exclusion:
            from lib.withExclusion import simulation, mul_sim_last, mul_sim_all
            self._sim = simulation
            self._mul_sim = mul_sim_last
            self._mul_sim_all = mul_sim_all
            self.options = self.parameters.as_tuple_without_cooperativity()

        elif self.parameters.q == 1 and not self.exclusion:
            from lib.noExclusion import simulation, mul_sim_last, mul_sim_all
            self._sim = simulation
            self._mul_sim = mul_sim_last
            self._mul_sim_all = mul_sim_all
            self.options = self.parameters.as_tuple_without_cooperativity()

        elif self.parameters.q != 1 and self.exclusion:
            from lib.Co_withExclusion import simulation, mul_sim_last, mul_sim_all
            self._sim = simulation
            self._mul_sim = mul_sim_last
            self._mul_sim_all = mul_sim_all
            self.options = self.parameters.as_tuple()

        else:
            from lib.Co_noExclusion import simulation, mul_sim_last, mul_sim_all
            self._sim = simulation
            self._mul_sim = mul_sim_last
            self._mul_sim_all = mul_sim_all
            self.options = self.parameters.as_tuple()

    def update_catalog(self, id: int):
        """
        Updates the catalog with simulation information.

        Args:
            id (int): The identifier for the simulation.
        """
        new = {}
        options = self.parameters.as_dict()
        for option in options.keys():
            new[option] = [options[option]]

        new["id"] = id
        new["cooperativity"] = int((self.parameters.q == 1))
        new["exclusion"] = self.exclusion
        master = pd.read_csv("./catalog.csv")
        new = pd.DataFrame(new)
        master = pd.concat([master, new])
        master.to_csv("./catalog.csv", index=False)

    def configure_id(self) -> int:
        """
        Configures the identifier for the simulation.

        Returns:
            int: The updated identifier.
        """
        with open(Paths.last_id_path, "rb") as f:
            x = pickle.load(f)

        with open(Paths.last_id_path, "wb") as f:
            pickle.dump(x + 1, f)

        return x + 1

    def run_write(self, num_of_sims: int):
        """
        Runs the simulation and writes the results to a file.

        Args:
            num_of_sims (int): Number of simulations to run.
        """
        id = self.configure_id()

        # Warm_up for NJIT
        data, times = self._mul_sim_all(1, *self.options)

        data, times = self._mul_sim_all(num_of_sims, *self.options)

        data_folder = os.path.abspath(Paths.data_folder)
        output_file_path = os.path.join(data_folder, str(id))
        with open(output_file_path, "wb") as f:
            pickle.dump({"data": data, "times": times}, f)

        self.update_catalog(id)
