import numpy as np
import numpy
import sys


sys.path.insert(1, '../NMR')
sys.path.insert(2, '../Misc')

import geometry
import CSA
import EFG
from Utils import exceptions
from Utils import overload



class nuclei:
    """
    A class to represent a nuclei.
    Contains information about coordinates, CSA and EFG
    Initialized with the isotope, and optionnaly geometry, CSA and EFG.
    
    ...
    
    Parameters
    ----------
    
    Attributes
    ----------
    isotope : str
        isotope
    geometry : class geometry.coordinates
        coordinates of the nuclei
    CSA : class CSA.CSA
        CSA tensor parameters
    EFG : class EFG.EFG
        EFG asymmetry
    
    Methods
    -------
    add_coordinates(system_frame=str, coordinates=numpy.ndarray)
    add_coordinates(system_frame=str, coordinates=list)
    add_coordinates(coordinates=numpy.ndarray)
    add_coordinates(coordinates=nlist)
        adds coordinates
    change_coordinates(system_frame=str, coordinates=numpy.ndarray)
    change_coordinates(system_frame=str, coordinates=list)
    change_coordinates(coordinates=numpy.ndarray)
    change_coordinates(coordinates=nlist)
        changes coordinates
    add_CSA(symmetry=str, amplitude=float)
    add_CSA(symmetry=str, amplitude=int)
    add_CSA(symmetry=str, amplitude=numpy.ndarray, orientation=numpy.ndarray)
    add_CSA(symmetry=str, amplitude=float, orientation=numpy.ndarray)
    add_CSA(symmetry=str, amplitude=int, orientation=numpy.ndarray)
        adds CSA
    change_CSA(symmetry=str, amplitude=float)
    change_CSA(symmetry=str, amplitude=int)
    change_CSA(symmetry=str, amplitude=numpy.ndarray, orientation=numpy.ndarray)
    change_CSA(symmetry=str, amplitude=float, orientation=numpy.ndarray)
    change_CSA(symmetry=str, amplitude=int, orientation=numpy.ndarray)
        changes CSA
    remove_CSA()
        resets CSA
    add_EFG(value=float)
    add_EFG(value=int)
        adds EFG
    change_EFG(value=float)
    change_EFG(value=int)
        changes EFG
    remove_EFG()
        resets EFG
    """
    
    def __init__(self, isotope, system_frame = None, coordinates = None, CSA_symmetry = None, CSA_amplitude = None, CSA_orientation = None, EFG_value = None):
        exceptions.isotope_check(isotope)
        
        self.isotope = isotope
        
        if isinstance(coordinates, (list, numpy.ndarray)):
            if isinstance(system_frame, str):
                self.geometry = geometry.coordinates(system_frame, coordinates)
            else:
                self.geometry = geometry.coordinates(coordinates)
        else:
            if coordinates == None:
                self.geometry = geometry.coordinates()
            else:
                raise SyntaxError('Provide a system_frame (Cartesian or Polar) and coordinates (numpy array) to initialize coordinates.')
            
        if isinstance(system_frame, str) and isinstance(CSA_amplitude, (int, float, numpy.ndarray)):
            if isinstance(CSA_orientation, numpy.ndarray): 
                self.CSA = CSA.CSA(system_frame, CSA_amplitude, CSA_orientation)
            elif CSA_orientation == None:
                self.CSA = CSA.CSA(system_frame, CSA_amplitude)
            else:
                raise SyntaxError('CSA_orientation must be a numpy array if not None')
        else:
            if CSA_symmetry == None or CSA_amplitude == None or CSA_orientation == None:
                self.CSA = CSA.CSA()
            else:
                raise SyntaxError('Provide a CSA_symmetry (isotropic, axial or asymmetric), CSA_amplitude and CSA_orientation (except for isotropic CSA).')
                
        self.EFG = EFG.EFG(isotope, EFG_value)
        
    
    @overload.overload_class_module((str, numpy.ndarray,))
    def add_coordinates(self, system_frame, coordinates):
        """
        adds coordinates.

        Parameters
        ----------
        system_frame : TYPE: str (optional)
            DESCRIPTION: system frame used (Cartesian or Polar). Default is Cartesian.
        coordinates : TYPE: numpy.array or list
            DESCRIPTION: coordinates of the nuclei in the chosen frame.

        Returns
        -------
        None.

        """
        
        self.geometry.add_coordinates(system_frame, coordinates)
        
    @overload.overload_class_module((str, list,))
    def add_coordinates(self, system_frame, coordinates):
        self.geometry.add_coordinates(system_frame, np.asarray(coordinates))
        
    @overload.overload_class_module((numpy.ndarray,))
    def add_coordinates(self, coordinates):
        self.geometry.add_coordinates('Cartesian', coordinates)
        
    @overload.overload_class_module((list,))
    def add_coordinates(self, coordinates):
        self.geometry.add_coordinates('Cartesian', np.asarray(coordinates))
        
        
    @overload.overload_class_module((str, numpy.ndarray,))
    def change_coordinates(self, system_frame, coordinates):
        """
        changes the coordinates

        Parameters
        ----------
        system_frame : TYPE: str (optional)
            DESCRIPTION: system frame used (Cartesian or Polar). Default is Cartesian.
        coordinates : TYPE: numpy.array or list
            DESCRIPTION: coordinates of the nuclei in the chosen frame.

        Returns
        -------
        None.

        """
        
        self.geometry.change_coordinates(system_frame, coordinates)
        
    @overload.overload_class_module((str, list,))
    def change_coordinates(self, system_frame, coordinates):
        self.geometry.change_coordinates(system_frame, numpy.asarray(coordinates))
        
    @overload.overload_class_module((numpy.ndarray,))
    def change_coordinates(self, coordinates):
        self.geometry.change_coordinates('Cartesian', coordinates)
        
    @overload.overload_class_module((list,))
    def change_coordinates(self, coordinates):
        self.geometry.change_coordinates('Cartesian', numpy.asarray(coordinates))
        
        
    @overload.overload_class_module((str, float,))
    def add_CSA(self, symmetry, amplitude):
        """
        adds CSA parameters

        Parameters
        ----------
        symmetry : TYPE: str
            DESCRIPTION: symmetry of the tensor (isotropic, axial or asymmetric).
        amplitude : TYPE: numpy.array, float or int
            DESCRIPTION: amplitude of the CSA tensor(s).
        orientation : TYPE: numpy.array (optional)
            DESCRIPTION: orientation of the CSA tensor

        Returns
        -------
        None.

        """
        
        self.CSA.add_CSA(symmetry, amplitude)
        
    @overload.overload_class_module((str, int,))
    def add_CSA(self, symmetry, amplitude):
        self.CSA.add_CSA(symmetry, float(amplitude))
        
    @overload.overload_class_module((str, numpy.ndarray, numpy.ndarray,))
    def add_CSA(self, symmetry, amplitude, orientation):
        self.CSA.add_CSA(symmetry, amplitude, orientation)
        
    @overload.overload_class_module((str, float, numpy.ndarray,))
    def add_CSA(self, symmetry, amplitude, orientation):
        self.CSA.add_CSA(symmetry, np.asarray([amplitude]), orientation)
        
    @overload.overload_class_module((str, int, numpy.ndarray,))
    def add_CSA(self, symmetry, amplitude, orientation):
        self.CSA.add_CSA(symmetry, np.asarray([float(amplitude)]), orientation)
        
    
    @overload.overload_class_module((str, float,))
    def change_CSA(self, symmetry, amplitude):
        """
        changes CSA parameters

        Parameters
        ----------
        symmetry : TYPE: str
            DESCRIPTION: symmetry of the tensor (isotropic, axial or asymmetric).
        amplitude : TYPE: numpy.array, float or int
            DESCRIPTION: amplitude of the CSA tensor(s).
        orientation : TYPE: numpy.array (optional)
            DESCRIPTION: orientation of the CSA tensor

        Returns
        -------
        None.

        """
        
        self.CSA.change_CSA(symmetry, amplitude)
        
    @overload.overload_class_module((str, int,))
    def change_CSA(self, symmetry, amplitude):
        self.CSA.change_CSA(symmetry, float(amplitude))
        
    @overload.overload_class_module((str, numpy.ndarray, numpy.ndarray,))
    def change_CSA(self, symmetry, amplitude, orientation):
        self.CSA.change_CSA(symmetry, amplitude, orientation)
        
    @overload.overload_class_module((str, float, numpy.ndarray,))
    def change_CSA(self, symmetry, amplitude, orientation):
        self.CSA.change_CSA(symmetry, np.asarray([amplitude]), orientation)
        
    @overload.overload_class_module((str, int, numpy.ndarray,))
    def change_CSA(self, symmetry, amplitude, orientation):
        self.CSA.change_CSA(symmetry, np.asarray([float(amplitude)]), orientation)
        
    
    def remove_CSA(self):
        """
        reset the CSA tensor to default

        Parameters
        ----------

        Returns
        -------
        None.
        
        """
        
        self.CSA.remove_CSA()
        
        
    @overload.overload_class_module((float,))
    def add_EFG(self, value):
        """
        adds EFG values

        Parameters
        ----------
        value : TYPE: int/float
            DESCRIPTION.: value of the EFG tensor

        Returns
        -------
        None.

        """
        
        self.EFG.add_EFG(value)
        
    @overload.overload_class_module((int,))
    def add_EFG(self, value):
        self.EFG.add_EFG(float(value))
        
        
    @overload.overload_class_module((float,))
    def change_EFG(self, value):
        """
        changes EFG values

        Parameters
        ----------
        value : TYPE: int/float
            DESCRIPTION.: value of the EFG tensor

        Returns
        -------
        None.

        """
        
        self.EFG.change_EFG(value)
        
    @overload.overload_class_module((int,))
    def change_EFG(self, value):
        self.EFG.change_EFG(float(value))
        
        
    def remove_EFG(self):
        """
        reset the EFG tensor to default

        Parameters
        ----------

        Returns
        -------
        None.
        
        """
        
        self.EFG.remove_EFG()