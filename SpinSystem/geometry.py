import numpy as np
import numpy
import sys

sys.path.insert(1, '../Misc')

from Utils import exceptions
from Utils import overload



class coordinates:
    """
    A class to represent the geometry of the spin system.
    Initialized with a dictionnary of nucleus: isotope.
    
    ...
    
    Parameters
    ----------
    system : TYPE: str, optional
        DESCRIPTION: system frame used (Cartesian or Polar).
    coordinates : TYPE: list/array, optional
        DESCRIPTION: coordinates of the nuclei.
    
    Attributes
    ----------
    coordinates : array
        coordinates of the nucleus
    
    Methods
    -------
    add_coordinates(system=str, Vector=numpyarray)
    add_coordinates(system=str, Vector=list)
    add_coordinates(Vector=numpy.array)
    add_coordinates(Vector=list)
        adds coordinates
    change_coordinates(system=str, Vector=numpy.ndarray)
    change_coordinates(system=str, Vector=list)
    change_coordinates(Vector=numpy.ndarray)
    change_coordinates(Vector=list)
        changes coordinates
    remove_coordinate()
        resets coordinates
    """
    
    @overload.overload_class_module(())
    def __init__(self):
        self.coordinates = []
        
    @overload.overload_class_module((str, numpy.ndarray,))
    def __init__(self, system, coordinates):
        self.coordinates = []
        self.add_coordinates(system, coordinates)
        
    @overload.overload_class_module((str, list))
    def __init__(self, system, coordinates):
        self.coordinates = []
        self.add_coordinates(system, coordinates)
        
        
    @overload.overload_class_module((str, numpy.ndarray,))
    def add_coordinates(self, system, Vector):
        """
        adds coordinates. 

        Parameters
        ----------
        system : TYPE: str, optional
            DESCRIPTION: axis systems used (Cartesian or Polar)
        Vector : TYPE: list/array or dict
            DESCRIPTION: coordinates in the choosen frame

        Raises
        ------
        ValueError
            DESCRIPTION: raised if coordinates are already defined

        Returns
        -------
        None.

        """
        
        if self.coordinates:
            raise ValueError('Coordinates defined. Use change_coordinate to modify them.')
        self.change_coordinates(system, Vector)
        
    @overload.overload_class_module((str, list,))
    def add_coordinates(self, system, Vector):
        self.add_coordinates(system, np.asarray(Vector))
        
    @overload.overload_class_module((numpy.ndarray,))
    def add_coordinates(self, Vector):
        self.add_coordinates('Cartesian', Vector)
        
    @overload.overload_class_module((list,))
    def add_coordinates(self, Vector):
        self.add_coordinates('Cartesian', np.asarray(Vector))
            

    @overload.overload_class_module((str, numpy.ndarray))
    def change_coordinates(self, system, Vector):
        """
        changes the coordinates.

        Parameters
        ----------
        system : TYPE: str
            DESCRIPTION: axis system used (Cartesian or Polar).
        Vector : TYPE: list/array of size 3, containing floats
            DESCRIPTION: coordinates in the choosen frame.

        Raises
        ------
        IndexError
            DESCRIPTION: raised if Vector does not contain 3 values.

        Returns
        -------
        None.

        """
        
        exceptions.member_of_list_check(system, ['Cartesian', 'Polar'])
        
        if len(Vector) != 3:
            if system == 'Cartesian':
                raise IndexError('{Vector} must contain 3 values: x, y and z in the Cartesian axis frame')
            if system == 'Polar':
                raise IndexError('{Vector} must contain 3 values: theta, phi and r in the Polar axis frame')
        for coord in Vector:
            exceptions.type_check(coord, float)
        
        if system == 'Cartesian':
            x, y, z = float(Vector[0]), float(Vector[1]), float(Vector[2])
 
        elif system == 'Polar':
            theta, phi, r = float(Vector[0]), float(Vector[1]), float(Vector[2])
            x, y, z = r * np.sin(theta) * np.cos(phi), r * np.sin(theta) * np.sin(phi), r * np.cos(theta)
        
        self.coordinates = np.asarray([x, y, z])
        
    @overload.overload_class_module((str, list))
    def change_coordinates(self, system, Vector):
        self.change_coordinates(system, np.asarray(Vector))
        
    @overload.overload_class_module((list,))
    def change_coordinates(self, Vector):
        self.change_coordinates('Cartesian', Vector)
        
    @overload.overload_class_module((numpy.ndarray,))
    def change_coordinates(self, Vector):
        self.change_coordinates('Cartesian', Vector)
                

    def remove_coordinate(self):
        """
        removes coordinates.

        Parameters
        ----------

        Returns
        -------
        None.

        """
        
        self.coordinates = []