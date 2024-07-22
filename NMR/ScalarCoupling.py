import numpy as np
import numpy
import sys

sys.path.insert(1, '../Misc')

from Utils import exceptions
from Utils import overload



class scalar_coupling:
    """
    Scalar coupling class. Initialized using a list of nuclei,
    and optionaly values of coupling constants.
    
    ...
    
    Parameters
    ----------
    nuclei : TYPE: aray
        DESCRIPTION: nuclei labels.
    coupling_matrix : TYPE: float or array, optional
        DESCRIPTION: coupling constants
        The default is None.
            
    Attributes
    ----------
    coupling_constants : numpy 2D array
        coupling constant matrix
    __nuclei : numpy array
        nuclei labels
    
    Methods
    -------
    add_coupling(coupling_matrix=float/array, nuclei_1 = None, nuclei_2 = None)
        adds coupling. Possible to do the whole spin system at one
    replace_coupling(coupling_constant=float, nuclei_1=str, nuclei_2=str)
        replaces the coupling constant between two nuclei
    remove_nuclei(nuclei = str)
        removes a nuclei from the class
    add_nuclei(nuclei = str)
        adds a nuclei from the class
    get_coupling_constant(nuclei_1=str, nuclei_2=str)
        returns the coupling constant between two nuclei
    """
    
    @overload.overload_class_module((numpy.ndarray,))
    def __init__(self, nuclei):
        for atom in nuclei:
            exceptions.type_check(atom, str)
        exceptions.duplicate_check(nuclei)
        
        self.__nuclei = np.asarray(nuclei)
        self.coupling_constants = np.zeros( (len(nuclei), len(nuclei)) )
        
    @overload.overload_class_module((list,))
    def __init__(self, nuclei):
        self.__init__(np.asarray(nuclei))
            
    @overload.overload_class_module((numpy.ndarray, float,))
    def __init__(self, nuclei, coupling_matrix):
        for atom in nuclei:
            exceptions.type_check(atom, str)
        exceptions.duplicate_check(nuclei)
        
        self.__nuclei = np.asarray(nuclei)
        self.coupling_constants = np.zeros( (len(nuclei), len(nuclei)) )
        
        self.add_coupling(coupling_matrix)
        
    @overload.overload_class_module((numpy.ndarray, int,))
    def __init__(self, nuclei, coupling_matrix):
        self.__init__(nuclei, float(coupling_matrix))
        
    @overload.overload_class_module((list, float,))
    def __init__(self, nuclei, coupling_matrix):
        self.__init__(np.asarray(nuclei), coupling_matrix)
        
    @overload.overload_class_module((list, int,))
    def __init__(self, nuclei, coupling_matrix):
        self.__init__(np.asarray(nuclei), float(coupling_matrix))
        
    @overload.overload_class_module((numpy.ndarray, numpy.ndarray,))
    def __init__(self, nuclei, coupling_matrix):
        for atom in nuclei:
            exceptions.type_check(atom, str)
        exceptions.duplicate_check(nuclei)
        
        self.__nuclei = np.asarray(nuclei)
        self.coupling_constants = np.zeros( (len(nuclei), len(nuclei)) )
        
        self.add_coupling(coupling_matrix)
        
    @overload.overload_class_module((numpy.ndarray, list,))
    def __init__(self, nuclei, coupling_matrix):
        self.__init__(nuclei, np.asarray(coupling_matrix))
        
    @overload.overload_class_module((list, numpy.ndarray,))
    def __init__(self, nuclei, coupling_matrix):
        self.__init__(np.asarray(nuclei), coupling_matrix)
        
    @overload.overload_class_module((list, list,))
    def __init__(self, nuclei, coupling_matrix):
        self.__init__(np.asarray(nuclei), np.asarray(coupling_matrix))
            
            
    def __get_position_nuclei(self, nuclei):
        """
        get the position of a nuclei in the coupling constant matrix

        Parameters
        ----------
        nuclei : TYPE: str
            DESCRIPTION: nuclei in the spin system.

        Returns
        -------
        position : TYPE: int
            DESCRIPTION: position of the nuclei  in the coupling constant matrix.

        """
        
        exceptions.nuclei_is_defined(self.__nuclei, nuclei)
        
        position = np.where(np.asarray(list(self.__nuclei)) == nuclei)[0][0]
        
        return position
        
        
    def add_coupling(self, coupling_matrix, nuclei_1 = None, nuclei_2 = None):
        """
        adds coupling constant values

        Parameters
        ----------
        coupling_matrix : TYPE: int, float or 2D numpy array
            DESCRIPTION: value(s) of coupling constants.
        nuclei_1 : TYPE: str, optional
            DESCRIPTION: nuclei 1 for which coupling constant is defined. The default is None.
        nuclei_2 : TYPE: str, optional
            DESCRIPTION: nuclei 2 for which coupling constant is defined. The default is None.

        Raises
        ------
        SyntaxError
            DESCRIPTION: raised if only one nuclei is given (instead of 0 or 2),
            if two nuclei are given but coupling_matrix is not an int/float,
            or if the dimensions of coupling_matrix do not match the number of nuclei.
        ValueError
            DESCRIPTION: raised if coupling_matrix is an array and nuclei are specified, or
            if a coupling constant between nuclei_1 and nuclei_2 is already defined.

        Returns
        -------
        None.

        """
        
        exceptions.type_check(coupling_matrix, (int, float, numpy.ndarray))
        
        if isinstance(coupling_matrix, (int, float)):
            if nuclei_1 != None and nuclei_2 == None:
                raise SyntaxError('Provide two nuclei to define a single coupling constant')
            if nuclei_2 != None and nuclei_1 == None:
                raise SyntaxError('Provide two nuclei to define a single coupling constant')
                
            if nuclei_1 == None:
                if not np.any(coupling_matrix):
                    for nuclei_1 in range(len(self.__nuclei) - 1):
                        for nuclei_2 in range(nuclei_1 + 1, len(self.__nuclei)):
                            self.change_coupling(coupling_matrix, self.__nuclei[nuclei_1], self.__nuclei[nuclei_2])
                else:
                    raise ValueError('Specify the nuclei for which the coupling constant is added.')
                    
            else:
                if self.get_coupling_constant(nuclei_1, nuclei_2) != 0.0:
                    raise ValueError(f'{nuclei_1} and {nuclei_2} already have a coupling constant defined. Use change_coupling to change it.')
                self.change_coupling(coupling_matrix, nuclei_1, nuclei_2)
            
        else:
            if len(coupling_matrix) != len(self.__nuclei):
                raise SyntaxError(f'Either provide a single coupling constant, or an 2D-array of dimension ({len(self.nuclei)}, {len(self.nuclei)})')
            if len(np.shape(coupling_matrix)) != 2 or np.shape(coupling_matrix)[0] != np.shape(coupling_matrix)[1]:
                raise SyntaxError(f'Either provide a single coupling constant, or an 2D-array of dimension ({len(self.nuclei)}, {len(self.nuclei)})')

            if nuclei_1 != None or nuclei_2 != None:
                raise ValueError('Conflicting instructions: provide either a single coupling constant for specific nuclei, or a 2D-array of dimension ({len(self.nuclei)}, {len(self.nuclei)})')  
                
            for nuclei_1 in range(len(self.__nuclei) - 1):
                for nuclei_2 in range(nuclei_1 + 1, len(self.__nuclei)):
                    self.change_coupling(coupling_matrix[nuclei_1, nuclei_2], self.__nuclei[nuclei_1], self.__nuclei[nuclei_2])
                    
                    
    def change_coupling(self, coupling_constant, nuclei_1, nuclei_2):
        """
        replaces the value of a coupling constant

        Parameters
        ----------
        coupling_constant : TYPE: int/float
            DESCRIPTION: coupling constant value.
        nuclei_1 : TYPE: str
            DESCRIPTION: nuclei in the spin system.
        nuclei_2 : TYPE: str
            DESCRIPTION: nuclei in the spin system.
            
        Raises
        ------
        SyntaxError
            DESCRIPTION: raised if nuclei_1 = nuclei_2.

        Returns
        -------
        None.

        """
        
        exceptions.type_check(coupling_constant, float)
        exceptions.nuclei_is_defined(self.__nuclei, nuclei_1)
        exceptions.nuclei_is_defined(self.__nuclei, nuclei_2)
        
        if nuclei_1 == nuclei_2:
            raise SyntaxError('Coupling constants can only be defined for different nuclei')
            
        position_nuclei_1 = self.__get_position_nuclei(nuclei_1)
        position_nuclei_2 = self.__get_position_nuclei(nuclei_2)
        
        self.coupling_constants[position_nuclei_1, position_nuclei_2] = coupling_constant
        self.coupling_constants[position_nuclei_2, position_nuclei_1] = coupling_constant
        
        
    def remove_nuclei(self, nuclei):
        """
        removes a nuclei from the coupling matrix

        Parameters
        ----------
        nuclei : TYPE: str
            DESCRIPTION: nuclei of the spin system.

        Returns
        -------
        None.

        """
        
        exceptions.nuclei_is_defined(self.__nuclei, nuclei)
        
        position = np.where(np.asarray(list(self.__nuclei)) == nuclei)[0][0]
        
        coupling_input = self.coupling_constants.copy()
        self.coupling_constants = np.zeros( (len(self.__nuclei) - 1, len(self.__nuclei) - 1) )
        self.coupling_constants[:position, :position] = coupling_input[:position, :position]
        self.coupling_constants[position:, position:] = coupling_input[position + 1:, position + 1:]
        
        self.__nuclei = np.asarray(list(filter(lambda x: x != nuclei, self.__nuclei)))
        
        
    def add_nuclei(self, nuclei):
        """
        adds a nuclei to the coupling matrix

        Parameters
        ----------
        nuclei : TYPE: str
            DESCRIPTION: nuclei label.

        Returns
        -------
        None.

        """
        exceptions.type_check(nuclei, str)
        
        self.__nuclei = np.append(self.__nuclei, nuclei)
        
        coupling_new = np.zeros( (len(self.__nuclei), len(self.__nuclei)) )
        coupling_new[:-1, :-1] = self.coupling_constants
        
        self.coupling_constants = coupling_new
        
        
    def get_coupling_constant(self, nuclei_1, nuclei_2):
        """
        get the coupling constant between two nuclei

        Parameters
        ----------
        nuclei_1 : TYPE: str
            DESCRIPTION: nuclei in the spin system.
        nuclei_2 : TYPE: str
            DESCRIPTION: nuclei in the spin system.

        Returns
        -------
        J : TYPE: float
            DESCRIPTION: coupling constant between nuclei_1 and nuclei_2.

        """
        
        exceptions.nuclei_is_defined(self.__nuclei, nuclei_1)
        exceptions.nuclei_is_defined(self.__nuclei, nuclei_2)
        
        position_nuclei_1 = self.__get_position_nuclei(nuclei_1)
        position_nuclei_2 = self.__get_position_nuclei(nuclei_2)
        
        J = self.coupling_constants[position_nuclei_1, position_nuclei_2]
        
        return J