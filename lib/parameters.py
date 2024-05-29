class Parameters:
    """
    Class representing parameters for a simulation.

    Attributes:
        L (int): Length of the lattice.
        T (int): Total time steps for the simulation.
        kon (float): Association rate constant.
        koff (float): Dissociation rate constant.
        kstep (float): Time step for the simulation.
        q (float): Surface concentration of receptors.
        kq (float): On-rate constant.
    """

    def __init__(self, L: int, T: int, kon: float, koff: float, kstep: float, q: float, kq: float):
        """
        Initialize parameters with given parameters.

        Args:
            L (int): Length of the lattice.
            T (int): Total time steps for the simulation.
            kon (float): Association rate constant.
            koff (float): Dissociation rate constant.
            kstep (float): Time step for the simulation.
            kq (float): On-rate constant.
            q (float): Surface concentration of receptors.
        """
        self.L = L
        self.T = T
        self.kon = kon
        self.koff = koff
        self.kstep = kstep
        self.kq = kq
        self.q = q

    def __str__(self):
        """
        Return a string representation of the parameters object with each attribute on a separate line.
        """
        return f"L: {self.L}\nT: {self.T}\nkon: {self.kon}\nkoff: \
            {self.koff}\nkstep: {self.kstep}\nkq: {self.kq}\nq: {self.q}"
    
    def as_dict(self):
        """
        Return a dictionary representation of the parameters object with variable names as keys.
        """
        return {attr: getattr(self, attr) for attr in vars(self)}
    
    def as_tuple(self):
        """
        Return a tuple of all the parameters (L, T, kon, koff, kstep, kq, q)
        """
        return tuple(getattr(self, attr) for attr in vars(self))
    

    def as_tuple_without_cooperativity(self):
        """
        Return a tuple of all the parameters (L, T, kon, koff, kstep, kq, q)
        """
        return tuple(getattr(self, attr) for attr in vars(self) if attr not in ['kq', 'q'])
    
    
    def are_valid_types(self)->bool:
        """
        Returns True if the given parameters have correct types for a simulation.
        """
        # Define the expected types for each attribute
        expected_types = {
            'L': int,
            'T': int,
            'kon': float,
            'koff': float,
            'kstep': float,
            'q': float,
            'kq': float
        }

        # Check if the types match the expected types
        for attr, expected_type in expected_types.items():
            if not isinstance(getattr(self, attr), expected_type):
                return False

        return True

    def mathematically_valid_ranges(self)->bool:
        """ 
        Returns True if the given options are in the 
        correct mathematical range.
        """

        for attr in vars(self):
            value = getattr(self, attr)
            if value < 0:
                return False
            if value == 0 and attr not in ['q', 'kq']:
                return False
        
        return True
    

    def is_q_and_kq_valid(self)->bool:
        """ 
        Return True if both q and kq are both 0 or non zero
        """
        return (self.q==1) == (self.kq==0)
    

    def no_debug_check(self)->bool:
        if not self.are_valid_types:
            return False
        if not self.mathematically_valid_ranges:
            return False
        if not self.is_q_and_kq_valid:
            return False
        
        return True
    
    def debug_check(self)->bool|None:
        if not self.are_valid_types:
            raise TypeError("The types of the parameters are not correct.")
        if not self.mathematically_valid_ranges:
            raise ValueError("The passed numbers are not mathematically correct.")
        if not self.is_q_and_kq_valid:
            raise ValueError("q and kq are not set correctly. if q=1 kq can only be 0")
        
        return True
    
    def check(self, debug:bool=False)->bool:
        if debug:
            return self.debug_check()
        else:
            return self.no_debug_check()