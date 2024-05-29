
import numpy as np
import pickle
from typing import List, Tuple
from lib.parameters import Parameters
import pandas as pd
from lib.paths import Paths
import os
from lib.simulation import Simulation, run_simulations


class Reader:

    """ 
    Reader class:
    The Reader class is used to load and read simulation data based on parameters.

Methods:
    __init__(parameters: Parameters)
        Initializes the Reader object with the provided Parameters object.
        If no matching simulation file is found for the given parameters, it raises a FileExistsError.

    __str__()
        Returns a string representation of the Reader object, including its parameters, ID, and path.

    set_id(id)
        Sets the ID of the Reader object and loads the associated path.

    load()
        Loads the data and times from the simulation file associated with the Reader object.
        If no matching simulation file is found, it raises a FileExistsError.

    get_last_frames()
        Returns the last frames of the simulation data.

    _load_pickle()
        Loads data and times from the simulation file.

    _load_path()
        Constructs the path of the simulation file based on the ID.

    _load_id()
        Loads the ID of the simulation based on the parameters.

    _get_matching_ids()
        Returns IDs of simulations matching the parameters.

    _do_par_match_id() -> bool
        Checks if the ID matches the parameters.

    _is_file_exist() -> bool
        Checks if a simulation file exists for the given parameters.
    
    """

    def __init__(self, parameters: Parameters):
        """
        Initializes the Reader object.

        Args:
            parameters (Parameters): An instance of the Parameters class.
                                     Contains parameters of the simulation.
        """
        self.parameters: Parameters = parameters
        self.id = None
        self.path = None
        self.data = None
        self.times = None

        if not self._is_file_exist():
            raise FileExistsError("Parameters do not match with a simulation. "
                                  "Run one first.")

    def __str__(self):
        """
        Returns a string representation of the Reader object.

        Returns:
            str: A string containing parameters, ID, and path of the Reader object.
        """
        return f" ----Parameters-------- \n {self.parameters}\n" \
               f" --------------------- \nID: {self.id}\nPATH: {self.path}\n"

    def set_id(self, id):
        """
        Sets the ID of the Reader object and loads associated path.

        Args:
            id: The ID of the simulation.
        """
        self.id = id
        self._load_path()
        if not self._do_par_match_id():
            raise ValueError("This value does not correspond to a simulation in catalog.")

    def _load(self):
        """
        Loads the data and times from the simulation file associated with the Reader object.

        Raises:
            FileNotFoundError: If no matching simulation file is found.
        """
        if not self._is_file_exist():
            raise FileNotFoundError("No simulation file found for the given parameters. Run a simulation first.")

        self._load_id()
        self._load_path()
        self._load_pickle()


    def load(self):
        """
        Loads the data and times from the simulation file associated with the Reader object.
        Assumes that you know what you are doing.
        
        Raises:
            FileNotFoundError: If no matching simulation file is found.
        """
        self._load_pickle()


    def get_last_frames(self):
        """
        Returns the last frames of the simulation data.

        Returns:
            numpy.ndarray: Last frames of the simulation data.
        """
        if self.data == None:
            self.load()
        return self.data[:, -1, :]

    def _load_pickle(self):
        """
        Loads data and times from the simulation file.

        Raises:
            ValueError: If the loaded data does not contain expected keys.
        """
        with open(self.path, "rb") as f:
            dic = pickle.load(f)
            if "data" not in dic or "times" not in dic:
                raise ValueError("Loaded data does not contain expected keys.")
            self.data = dic["data"]
            self.times = dic["times"]


    def _load_path(self):
        """
        Constructs the path of the simulation file based on the ID.
        """
        self.path = os.path.join(Paths.data_folder, str(self.id))

    def _load_id(self):
        """
        Loads the ID of the simulation based on the parameters.
        """
        filtered_df = pd.read_csv(Paths.catalog)
        my_dict = self.parameters.as_dict()

        for key, value in my_dict.items():
            filtered_df = filtered_df[filtered_df[key] == value]

        self.id = filtered_df["id"].values[0]

    def get_matching_ids(self):
        """
        Returns IDs of simulations matching the parameters.

        Returns:
            numpy.ndarray: Array of matching IDs.
        """
        filtered_df = pd.read_csv(Paths.catalog)
        my_dict = self.parameters.as_dict()

        for key, value in my_dict.items():
            filtered_df = filtered_df[filtered_df[key] == value]

        return filtered_df["id"].values
    
    def _get_filtered_dataframe(self) -> pd.DataFrame:
        """
        Returns a filtered DataFrame based on the parameters.

        Returns:
            pd.DataFrame: Filtered DataFrame containing matching simulation entries.

        Raises:
            FileNotFoundError: If the catalog file does not exist.
        """
        if not os.path.exists(Paths.catalog):
            raise FileNotFoundError("Catalog file not found.")

        if not hasattr(self, '_filtered_df'):
            filtered_df = pd.read_csv(Paths.catalog)
            my_dict = self.parameters.as_dict()

            for key, value in my_dict.items():
                filtered_df = filtered_df[filtered_df[key] == value]
            
            self._filtered_df = filtered_df

        return self._filtered_df



    def _do_par_match_id(self) -> bool:
        """
        Checks if the ID matches parameters.

        Returns:
            bool: True if ID matches parameters, False otherwise.
        """
        return (self.id in self.get_matching_ids())

    def _is_file_exist(self) -> bool:
        """
        Checks if a simulation file exists for the given parameters.

        Returns:
            bool: True if a simulation file exists, False otherwise.
        """
        filtered_df = pd.read_csv(Paths.catalog)
        my_dict = self.parameters.as_dict()

        for key, value in my_dict.items():
            filtered_df = filtered_df[filtered_df[key] == value]

        if len(filtered_df) == 1:
            return True
        elif len(filtered_df) > 1:
            print("Multiple simulation files found. Please set the ID explicitly using set_id()\n \
                    Use get_matching_ids() to get all the ids that are compatible with these parameters.")
            return True
        return False