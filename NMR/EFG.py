import sys

sys.path.insert(1, '../Misc')

import constants as cst
from Utils import exceptions
from Utils import overload



class EFG:
    """
    A class to represent the Electric Field Gradient.
    Initialized with an isotope, and optionally specifiying EFG parameters
    
    ...
    
    Parameters
    ----------
    isotope : TYPE: str
        DESCRIPTION: isotope.
    EFG_input : TYPE: int or float, optional
        DESCRIPTION: EFG value.
        The default is None.
    
    Attributes
    ----------
    EFG_asymmetry : float
        asymetry of the EFG. Default is None
    __isotope : str
        isotope
    
    Methods
    -------
    add_EFG(value=float)
    add_EFG(value=int)
        adds a EFG value
    change_EFG(value=float)
    change_EFG(value=int)
        changes a EFG value
    remove_EFG()
        resets the EFG
    """
    
    def __init__(self, isotope, EFG_input = None):
        exceptions.isotope_check(isotope)
        
        if EFG_input == None:
            self.EFG_asymmetry = None
            self.__isotope = isotope
        else:
            exceptions.type_check(EFG_input, (int, float))
            self.__isotope = isotope
            
            if cst.spin_quantum_number[isotope] == '1/2':
                self.EFG_asymmetry = None
            else:
                self.EFG_asymmetry = float(EFG_input)

               
    @overload.overload_class_module((float,))
    def add_EFG(self, value):
        """
        sets the value of the EFG asymmetry of a nucleus

        Parameters
        ----------
        value : TYPE: int/float
            DESCRIPTION: value of the EFG.
            
        Raises
        ------
        ValueError
            DESCRIPTION: raised if nuclei is already has an EFG

        Returns
        -------
        None.

        """
        
        if self.EFG_asymmetry != None:
            raise ValueError('EFG already defined. Use change_EFG to modify the EFG tensor.')
            
        self.change_EFG(value)
        
    @overload.overload_class_module((int,))
    def add_EFG(self, value):
        self.add_EFG(float(value))
        
        
    @overload.overload_class_module((float,))
    def change_EFG(self, value):
        """
        changes the value of the EFG asymmetry of a nucleus

        Parameters
        ----------
        value : TYPE: int/float
            DESCRIPTION: value of the EFG.

        Raises
        ------
        SyntaxError
            DESCRIPTION: raised if trying to set the EFG for a spin-1/2 nucleus.

        Returns
        -------
        None.

        """
        
        if cst.spin_quantum_number[self.__isotope] == '1/2':
            raise SyntaxError('Asymmetry of EFG tensor can only be set for non spin-1/2 nuclei')
        
        self.EFG_asymmetry = value
        
    @overload.overload_class_module((int,))
    def change_EFG(self, value):
        self.change_EFG(float(value))
        
        
    def remove_EFG(self):
        """
        resets EFG

        Parameters
        ----------

        Returns
        -------
        None.

        """
        
        self.EFG_asymmetry = None