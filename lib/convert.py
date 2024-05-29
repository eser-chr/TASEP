class SideLengthConverter:
    """
    Utility class for converting between side length and distance units.
    """
    _SIDE_LENGTH: float = 8e-3  # Default side length in meters (8nm)

    @staticmethod
    def set_side_length(side_length: float) -> None:
        """
        Set the side length for conversion.

        Args:
            side_length (float): The side length in meters.
        """
        SideLengthConverter._SIDE_LENGTH = side_length

    @staticmethod
    def convert_sides_to_distance_mu(num_sides: int) -> float:
        """
        Convert number of sides to distance in micrometers.

        Args:
            num_sides (int): Number of sides.

        Returns:
            float: Distance in micrometers.
        """
        return SideLengthConverter._SIDE_LENGTH * num_sides

    @staticmethod
    def convert_sides_to_distance_nm(num_sides: int) -> float:
        """
        Convert number of sides to distance in nanometers.

        Args:
            num_sides (int): Number of sides.

        Returns:
            float: Distance in nanometers.
        """
        return SideLengthConverter._SIDE_LENGTH * num_sides * 1000

    @staticmethod
    def convert_distance_mu_to_sides(length: float) -> int:
        """
        Convert distance in micrometers to number of sides.

        Args:
            length (float): Distance in micrometers.

        Returns:
            int: Number of sides.
        """
        return int(length / SideLengthConverter._SIDE_LENGTH)

    @staticmethod
    def convert_distance_nm_to_sides(length: float) -> int:
        """
        Convert distance in nanometers to number of sides.

        Args:
            length (float): Distance in nanometers.

        Returns:
            int: Number of sides.
        """
        return int(length / (1000 * SideLengthConverter._SIDE_LENGTH))
