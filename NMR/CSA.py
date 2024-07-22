import numpy as np
import numpy
import sys

sys.path.insert(1, '../Misc')

from Utils import exceptions
from Utils import overload



def get_CSA_amplitude(amplitude_input, symmetry):
    """
    get the CSA amplitude of CSA component(s).

    Parameters
    ----------
    amplitude_input : TYPE: float/int or list/array
        DESCRIPTION: amplitude.
    symmetry : TYPE: str
        DESCRIPTION: symmetry of the CSA tensory (isotopic, axial or asymmetric).

    Raises
    ------
    TypeError
        DESCRIPTION: raised if the amplitude_input is not a number/array.

    Returns
    -------
    csa_val : TYPE: float or list
        DESCRIPTION: amplitude(s) of the CSA component(s).

    """
    
    if exceptions.member_of_list_check(symmetry, ['isotopic', 'axial', 'asymmetric']):
        if symmetry == 'isotropic' or symmetry == 'axial':
            if isinstance(amplitude_input, (int, float)):
                csa_val = float(amplitude_input)
            elif isinstance(amplitude_input, (list, numpy.ndarray)) and np.shape(np.asarray(amplitude_input)) == (1,):
                csa_val = float(amplitude_input[0])
            else:
                raise TypeError(f'{amplitude_input} must be set to a floating number. \
                                CSA not changed.')
            
        elif symmetry == 'asymmetric':
            if isinstance(amplitude_input, (list, numpy.ndarray)) and len(amplitude_input) == 2:
                for amp in amplitude_input:
                    exceptions.type_check(amp, float)
            else:
                raise TypeError('An array of CSA values must be provided. \
                                  CSA not changed.')
            csa_val = [amplitude_input[0], amplitude_input[1]]
        
    return csa_val


def get_CSA_component_orientation(orientation_input):
    """
    get the CSA orientation of a single CSA component in a cartesian coordinate.

    Parameters
    ----------
    orientation_input : TYPE: list or array of size 2 or 3
        DESCRIPTION: contains the orientation, either polar or cartesian.

    Raises
    ------
    IndexError
        DESCRIPTION: raised if size of orientaiton_input is not 2 or 3.

    Returns
    -------
    orientation : TYPE: numpy array
        DESCRIPTION: cartesian orientation.

    """
    if len(orientation_input) == 3:
        for coord in orientation_input:
            exceptions.type_check(coord, float)
            
        orientation = np.asarray(orientation_input)
        return orientation
    
    if len(orientation_input) == 2:
        for angle in orientation_input:
            exceptions.type_check(angle, float)
            
        theta, phi = orientation_input
        x, y, z = np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)
        orientation = np.asarray([x, y, z])
        return orientation
        
    raise IndexError(f'{orientation_input} should be of size 2 (polar angle) or 3 (cartesian vector coordinate).\
                      CSA not changed.')
    

def get_CSA_orientation(orientation_input, symmetry):
    """
    gets the CSA orientation for all type of CSA tensor symmetry.

    Parameters
    ----------
    orientation_input : TYPE: list/array
        DESCRIPTION: contains the orientation of 1 or 2 CSA component.
    symmetry : TYPE: str
        DESCRIPTION: symmetry of the CSA tensor.

    Raises
    ------
    IndexError
        DESCRIPTION: raised when symmetry=asymmetric and len(orientation_input)!=2.
    TypeError
        DESCRIPTION: raised when orientation_input is not an array or list.
    ValueError
        DESCRIPTION; raised when symmetry is not axial or asymmetric.

    Returns
    -------
    orientation : TYPE: list
        DESCRIPTION: orientation of the CSA component(s).

    """
    
    if symmetry == 'axial':
        orientation = get_CSA_component_orientation(orientation_input)
        return orientation
                         
    elif symmetry == 'asymmetric':
        if isinstance(orientation_input, (list, numpy.ndarray)):
            if len(orientation_input) != 2:
                raise IndexError('Two CSA components must be provided.')
        else:
            raise TypeError('An array of CSA orientations must be provided.')
                             
        orientation = [get_CSA_component_orientation(orientation_input[0]),
                        get_CSA_component_orientation(orientation_input[1])]
        return orientation
                         
    else:
        raise ValueError('symmetry must be set to axial or asymmetric')


def change_CSA_isotropic(amplitude=0.0):
    """
    changes CSA for a symmetric chemical shift tensor.

    Parameters
    ----------
    amplitude : TYPE: float or list/array of shape (1,). Default is None
        DESCRIPTION: amplitude of the CSA. Default is 0
        
    Returns
    -------
    csa : TYPE: dictionnary
        DESCRIPTION: csa parameters for an isotropic tensor.

    """
    
    csa_amp = get_CSA_amplitude(amplitude, 'isotropic')
    
    csa = {'symmetry': 'isotropic',
            'amplitude': csa_amp,
            'orientation': None}
    
    return csa
    
                    
def change_CSA_axial(amplitude, orientation):
    """
    changes the CSA for an axially symmetric CSA tensor.

    Parameters
    ----------
    amplitude : TYPE: float or list/array of shape (1,)
        DESCRIPTION: amplitude of the CSA.
    orientation : TYPE: list/array of size 2 or 3
        DESCRIPTION: orientation (polar or cartesian) of the CSA tensor.

    Returns
    -------
    csa : TYPE: dictionnary
        DESCRIPTION: csa parameters for an axially symmetric tensor.
        if the amplitude is 0, returns a isotropic CSA

    """
    
    csa_amp = get_CSA_amplitude(amplitude, 'axial')
    csa_orient = get_CSA_orientation(orientation, 'axial')
                        
    csa = {'symmetry': 'axial',
            'amplitude': csa_amp,
            'orientation': csa_orient}
    
    return csa


def change_CSA_asymmetric(amplitude, orientation):
    """
    changes the CSA for an fully asymmetric CSA tensor.

    Parameters
    ----------
    amplitude : TYPE: list/array of size 2
        DESCRIPTION: contains the amplitudes.
    orientation : TYPE: list/array of size 2
        DESCRIPTION: contains the orientations..

    Raises
    ------
    IndexError
        DESCRIPTION: raised if amplitudes or orientation do not contain exactly 2 elements.
    TypeError
        DESCRIPTION: raised if amplitude or orientaion are not list/array.

    Returns
    -------
    csa : TYPE: dictionnary
        DESCRIPTION: csa parameters for an fully asymmetric tensor.
        If one of the amplitude is 0, returns an axially symmetric tensor.

    """
                         
    sigma_1, sigma_2 = get_CSA_amplitude(amplitude, 'asymmetric')
    orient_1, orient_2 = get_CSA_orientation(orientation, 'asymmetric')
    
    if sigma_1 == 0.0:
        return change_CSA_axial(sigma_2, orient_2)
    if sigma_2 == 0.0:
        return change_CSA_axial(sigma_1, orient_1)
    
    csa = {'symmetry': 'asymmetric',
            'amplitude': [sigma_1, sigma_2],
            'orientation': [orient_1, orient_2]}
    
    return csa



class CSA:
    """
    A class to represent the CSA tensors.
    Initialized by optionally specifiying CSA tensor parameters
    
    ...
    
    Parameters
    ----------
    symmetry : TYPE: str, optional
        DESCRIPTION: used for a deferent initialization of the symmetry.
    amplitude : TYPE: float or array, optional
        DESCRIPTION: used for a diferent initialization of the amplitudes.
    orientation : TYPE: array, optional
        DESCRIPTION: used for a diferent initialization of the orientations.
        
    Attributes
    ----------
    symmetry : str
        symmetry of the CSA tensor. Default is isotropic
    csa : float or array
        amplitude of the CSA tensor. Default is None
    orientation : array
        orientation of the CSA tensor. Default is None
    
    Methods
    -------
    add_CSA(ssymmetry=str, amplitude=float)
    add_CSA(ssymmetry=str, amplitude=int)
    add_CSA(ssymmetry=str, amplitude=float, orientation=numpy.ndarray)
    add_CSA(ssymmetry=str, amplitude=float, orientation=list)
    add_CSA(ssymmetry=str, amplitude=int, orientation=numpy.ndarray)
    add_CSA(ssymmetry=str, amplitude=int, orientation=list)
    add_CSA(ssymmetry=str, amplitude=numpy.ndarray, orientation=numpy.ndarray)
    add_CSA(ssymmetry=str, amplitude=list, orientation=numpy.ndarray)
    add_CSA(ssymmetry=str, amplitude=numpy.ndarray, orientation=nlist)
    add_CSA(ssymmetry=str, amplitude=list, orientation=list)
        adds CSA parameters
    change_CSA(symmetry=str, amplitude=float)
    change_CSA(symmetry=str, amplitude=int)
    change_CSA(symmetry=str, amplitude=float, orientation=numpy.ndarray)
    change_CSA(symmetry=str, amplitude=float, orientation=list)
    change_CSA(symmetry=str, amplitude=int, orientation=numpy.ndarray)
    change_CSA(symmetry=str, amplitude=int, orientation=list
    change_CSA(symmetry=str, amplitude=numpy.ndarray, orientation=numpy.ndarray)
    change_CSA(symmetry=str, amplitude=list, orientation=numpy.ndarray)
    change_CSA(symmetry=str, amplitude=numpy.ndarray, orientation=list)
    change_CSA(symmetry=str, amplitude=list, orientation=list)
        changes the CSA tensor parameters
    remove_CSA()
        resets the CSA parameters
    """
    
    @overload.overload_class_module(())
    def __init__(self):
        self.symmetry = 'isotropic'
        self.csa = None
        self.orientation = None
            
    @overload.overload_class_module((str, float))
    def __init__(self, symmetry, amplitude):
        self.symmetry = 'isotropic'
        self.csa = None
        self.orientation = None
        
        if symmetry != 'isotropic':
            print('Orientations must be provided. Default initialization.')
            
        self.symmetry = 'isotropic'
        self.csa = amplitude
        self.orientation = None
        
    @overload.overload_class_module((str, int))
    def __init__(self, symmetry, amplitude):
        self.__init__(symmetry, float(amplitude))
        
    @overload.overload_class_module((str, numpy.ndarray, numpy.ndarray,))
    def __init__(self, symmetry, amplitude, orientation):
        self.symmetry = 'isotropic'
        self.csa = None
        self.orientation = None
        
        if symmetry == 'isotropic':
            print('Orientations not needed for isotropic CSA. Default initialization.')
        exceptions.member_of_list_check(symmetry, ['isotropic', 'axial', 'asymmetric'])
            
        self.change_CSA(symmetry, amplitude, orientation)
            
    @overload.overload_class_module((str, float, numpy.ndarray,))
    def __init__(self, symmetry, amplitude, orientation):
        self.__init__(symmetry, np.asarray([amplitude]), orientation)
        
    @overload.overload_class_module((str, int, numpy.ndarray,))
    def __init__(self, symmetry, amplitude, orientation):
        self.__init__(symmetry, float(amplitude), orientation)
            
        
    @overload.overload_class_module((str, float))
    def add_CSA(self, symmetry, amplitude):
        """
        Adds a CSA value for a nucleus.
        This function uses the module change_CSA. Check associated docstring
        for more info.

        Parameters
        ----------
        symmetry : TYPE: str
            DESCRIPTION: symmetry of the CSA tensor.
        amplitude : TYPE: float or array
            DESCRIPTION: amplitude(s) of the CSA component.
        orientation : TYPE: array, optional
            DESCRIPTION: orientation(s) of the CSA tensor component (polar angle
                        of cartesian vector coordinate).
            
        Raises
        ------
        ValueError
            DESCRIPTION: raised if CSA is already defined

        Returns
        -------
        None.

        """
        
        if self.csa != None:
            raise ValueError('CSA defined. Use change_CSA to modify the tensor.')
            
        self.change_CSA(symmetry, amplitude)
        
    @overload.overload_class_module((str, int))
    def add_CSA(self, symmetry, amplitude):
        self.add_CSA(symmetry, float(amplitude))
        
    @overload.overload_class_module((str, numpy.ndarray, numpy.ndarray,))
    def add_CSA(self, symmetry, amplitude, orientation):
        if self.csa != None:
            raise ValueError('CSA defined. Use change_CSA to modify the tensor.')
            
        self.change_CSA(symmetry, amplitude, orientation)
        
    @overload.overload_class_module((str, list, numpy.ndarray,))
    def add_CSA(self, symmetry, amplitude, orientation):
        self.add_CSA(symmetry, np.asarray(amplitude), orientation)
        
    @overload.overload_class_module((str, numpy.ndarray, list,))
    def add_CSA(self, symmetry, amplitude, orientation):
        self.add_CSA(symmetry, amplitude, np.asarray(orientation))
        
    @overload.overload_class_module((str, list, list,))
    def add_CSA(self, symmetry, amplitude, orientation):
        self.add_CSA(symmetry, np.asarray(amplitude), np.asarray(orientation))
        
    @overload.overload_class_module((str, float, numpy.ndarray,))
    def add_CSA(self, symmetry, amplitude, orientation):
        self.add_CSA(symmetry, np.asarray([amplitude]), orientation)
        
    @overload.overload_class_module((str, float, list,))
    def add_CSA(self, symmetry, amplitude, orientation):
        self.add_CSA(symmetry, np.asarray([amplitude]), np.asarray(orientation))
        
    @overload.overload_class_module((str, int, numpy.ndarray,))
    def add_CSA(self, symmetry, amplitude, orientation):
        self.add_CSA(symmetry, float(amplitude), orientation)
        
    @overload.overload_class_module((str, int, list,))
    def add_CSA(self, symmetry, amplitude, orientation):
        self.add_CSA(symmetry, float(amplitude), np.asarray(orientation))

        
    @overload.overload_class_module((str, float))
    def change_CSA(self, symmetry, amplitude):
        """
        changes the characterics of the CSA tensor

        Parameters
        ----------
        symmetry : TYPE: str
            DESCRIPTION: symmetry of the CSA tensor.
        amplitude : TYPE: float or array
            DESCRIPTION: amplitude(s) of the CSA component.
        orientation : TYPE: array/list, optional
            DESCRIPTION: orientation(s) of the CSA tensor component (polar angle
                        of cartesian vector coordinate).

        Raises
        ------
        SyntaxError
            DESCRIPTION: raised if orientation is not provided and CSA is not isotropic.

        Returns
        -------
        None.

        """
        
        if symmetry != 'isotropic':
            raise SyntaxError('Orientations must be provided. Default initialization.')
        
        tensor = change_CSA_isotropic(amplitude)
                    
        self.symmetry = tensor['symmetry']
        self.csa = tensor['amplitude']
        self.orientation = None
        
    @overload.overload_class_module((str, int))
    def change_CSA(self, symmetry, amplitude):
        self.change_CSA(symmetry, float(amplitude))
        
    @overload.overload_class_module((str, numpy.ndarray, numpy.ndarray,))
    def change_CSA(self, symmetry, amplitude, orientation):
        if symmetry == 'isotropic':
            raise SyntaxError('Orientations not needed for isotropic CSA. CSA not changed.')
        exceptions.member_of_list_check(symmetry, ['isotropic', 'axial', 'asymmetric'])
        
        if symmetry == 'axial':
            tensor = change_CSA_axial(amplitude, orientation)
                
        elif symmetry == 'asymmetric':
            tensor = change_CSA_asymmetric(amplitude, orientation)
            
        self.symmetry = tensor['symmetry']
        self.csa = tensor['amplitude']
        self.orientation = tensor['orientation']
        
    @overload.overload_class_module((str, list, numpy.ndarray,))
    def change_CSA(self, symmetry, amplitude, orientation):
        self.change_CSA(symmetry, np.asarray(amplitude), orientation)
        
    @overload.overload_class_module((str, numpy.ndarray, list,))
    def change_CSA(self, symmetry, amplitude, orientation):
        self.change_CSA(symmetry, amplitude, np.asarray(orientation))
        
    @overload.overload_class_module((str, list, list,))
    def change_CSA(self, symmetry, amplitude, orientation):
        self.change_CSA(symmetry, np.asarray(amplitude), np.asarray(orientation))
        
    @overload.overload_class_module((str, float, numpy.ndarray,))
    def change_CSA(self, symmetry, amplitude, orientation):
        self.change_CSA(symmetry, np.asarray([amplitude]), orientation)
        
    @overload.overload_class_module((str, float, list,))
    def change_CSA(self, symmetry, amplitude, orientation):
        self.change_CSA(symmetry, np.asarray([amplitude]), np.asarray(orientation))
        
    @overload.overload_class_module((str, int, numpy.ndarray,))
    def change_CSA(self, symmetry, amplitude, orientation):
        self.change_CSA(symmetry, float(amplitude), orientation)
        
    @overload.overload_class_module((str, int, list,))
    def change_CSA(self, symmetry, amplitude, orientation):
        self.change_CSA(symmetry, float(amplitude), np.asarray(orientation))
        
        
    def remove_CSA(self):
        """
        reset the CSA tensor to default

        Parameters
        ----------

        Returns
        -------
        None.

        """
        
        self.symmetry = 'isotropic'
        self.csa = None
        self.orientation = None