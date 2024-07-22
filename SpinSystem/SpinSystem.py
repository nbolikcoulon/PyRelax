import numpy as np
import os
import sys

sys.path.insert(1, '../NMR')
sys.path.insert(2, '../Quantum')
sys.path.insert(3, '../Outputs')
sys.path.insert(4, '../Misc')

import nuclei
import Basis
import ScalarCoupling
import Text
import Figures
from Utils import exceptions
from Utils import overload



class spin_system:
    """
    A class to represent the spin system.
    Initialized with a dictionnary associating a label to an isotope for each
    nuclei in the spin system, in the form {label: isotope}, or {label: nuclei}
    
    ...
    
    Parameters
    ----------
    nuclei : TYPE: dictionnary
        DESCRIPTION: associate a label to a nuclei (class) for each nuclei in the
        spin system, in the form {label: nuclei}.
    
    Attributes
    ----------
    nuclei : dict of class nuclei.nuclei
        associates a nuclei to its coordinates, CSA and EFG parameters
    coupling_constants : class for the coupling constant matrix
        Default is coupling constants set to 0
    basis : class
        basis set of the spin system. Default is CartesianOperatorBasis
    equivalent_nuclei : dictionnary
        assocatiates a group to a list of equivalent nuclei. This is used 
        to symmetrize the basis when not empty
        Default is None
    name : str
        name of the spin system
    dir_out_name : str
        directory where outputs are stored
    
    Methods
    -------
    remove_nuclei(label=str):
        removes a nuclei from the spin system
    add_nuclei(label=str, isotope=str)
    add_nuclei(label=str, nuclei_input=nuclei.nuclei):
        adds a nuclei to the spin system
    spin_permutation(permutation=list/array)
        defines equivalent nuclei
    reinitialize_permutation()
        resets the spin symmetry
    show_spin_system(save=bool, figname=str, overwrite=bool)
        creates a figure of the geometry of the spin system
    Private methods
    ---------------
    __name_spin_system()
        name of the spin system
    __make_out_dir(name=str)
        creates the output directory
    """
    
    def __init__(self, nuclei_dict):
        exceptions.type_check(nuclei_dict, dict)
        
        type_input = list(nuclei_dict.values())[0]
        exceptions.type_check(type_input, (str, nuclei.nuclei))
        if type_input == str:
            exceptions.isotope_check(list(nuclei_dict.values()))
            
        self.nuclei = dict()
        for k,v in nuclei_dict.items():
            if isinstance(v, type(type_input)):
                self.nuclei[k] = nuclei.nuclei(v) if type(type_input) == str else v
            else:
                self.nuclei = None
                raise SyntaxError('nuclei_dict values should all have the same type (str or nuclei.nuclei).')
            
        self.coupling_constants = ScalarCoupling.scalar_coupling(list(self.nuclei.keys()))
        self.equivalent_nuclei = None
        self.basis = Basis.Basis(self.nuclei, 'CartesianOperatorBasis')
            
        self.name = self.__name_spin_system(self.nuclei)
        self.dir_out_name = self.__make_out_dir(self.name)
        print('Spin sytem set')
        print(f'Outputs will be placed in directory {self.dir_out_name[:-1]}')
        
        Text.write_spinsystem(self)
        
        
    @staticmethod
    def __name_spin_system(nuclei):
        """
        name of the spin system, in the form LaMbNc... with
        L, M, N nuclei of the spin system and a, b, c their multiplicity.
    
        Parameters
        ----------
        nuclei : TYPE: dict
            DESCRIPTION: associate a label to an isotope for each nuclei in the
            spin system, in the form {label: isotope}.
    
        Returns
        -------
        name : TYPE: str
            DESCRIPTION: name of the spin system, in the form LaMbNc... with
            L, M, N nuclei of the spin system and a, b, c their multiplicity.
    
        """
        
        counts = {}
        for isotope in nuclei.values():
            if isotope not in counts.keys():
                counts[isotope] = 0
            counts[isotope] += 1
            
        name = ''
        for isotope, count in counts.items():
            if count == 1:
                name += f'{isotope}_'
            else:
                name += f'{isotope}{count}_'
        name = name[:-1]
        
        return name
    
    
    @staticmethod
    def __make_out_dir(name):
        """
        creates output directory
    
        Parameters
        ----------
        name : TYPE: str
            DESCRIPTION: name of the spin system.
    
        Returns
        -------
        OutDir : TYPE: str
            DESCRIPTION: name of the ouput directory.
    
        """
        
        OutDir = os.getcwd() + '/' + name
        
        if not os.path.isdir(OutDir):
            os.mkdir(OutDir)
            return OutDir
        
        n = 1
        while os.path.isdir(OutDir):
            OutDir = name + '_' + str(n)
            n += 1
        os.mkdir(OutDir)
        return OutDir
    
            
    def remove_nuclei(self, label):
        """
        Removes a nuclei from the spin system.

        Parameters
        ----------
        label : TYPE: str
            DESCRIPTION: a nuclei in the spin system.

        Returns
        -------
        None.

        """
        
        exceptions.nuclei_is_defined(self.nuclei, label)
            
        if len(self.nuclei) == 1:
            raise IndexError('Cannot remove last nucleus from the spin system')

        del self.nuclei[label]
        
        self.coupling_constants.remove_nuclei(label)
        self.basis = Basis.Basis(self.nuclei, self.basis.basis_type)
        
        if self.equivalent_nuclei != None:
            for group, nuclei_list in self.equivalent_nuclei.items():
                if label in nuclei_list:
                    if len(nuclei_list) == 2:
                        del self.equivalent_nuclei[group]
                        if len(self.equivalent_nuclei) == 0:
                            self.equivalent_nuclei = None
                        break
                    else:
                        self.equivalent_nuclei[group] = np.asarray(list(filter(lambda x: x!=label, nuclei_list)))
                        break
        
        self.name = self.__name_spin_system(self.nuclei)
        self.dir_out_name = self.__make_out_dir(self.name)
        print('New spin sytem set')
        print(f'Outputs will be placed in directory {self.dir_out_name}')
        
        Text.write_spinsystem(self)

            
    @overload.overload_class_module((str, str))
    def add_nuclei(self, label, isotope):
        """
        adds a nuclei in the spin system.

        Parameters
        ----------
        Either:
            label : TYPE: str
                DESCRIPTION: atom label.
            isotope : TYPE: str
                DESCRIPTION: isotope.
        Or:
            label : TYPE: str
                DESCRIPTION: atom label.
            nuclei_input : TYPE: nuclei.nuclei
                DESCRIPTION: nuclei.
                
        Raises
        ------
        ValueError
            DESCRIPTION: raised if label is already defined in the spin system.

        Returns
        -------
        None.

        """
        
        if label in self.nuclei.keys():
            raise ValueError(f'{label} is already present in the spin system.')
        exceptions.isotope_check(isotope)
            
        self.nuclei[label] = nuclei.nuclei(isotope)
        self.coupling_constants.add_nuclei(label)
        self.basis = Basis.Basis(self.nuclei, self.basis.basis_type)
        
        self.name = self.__name_spin_system(self.nuclei)
        self.dir_out_name = self.__make_out_dir(self.name)
        print('New spin sytem set')
        print(f'Outputs will be placed in directory {self.dir_out_name}')
        
        Text.write_spinsystem(self)
        
    @overload.overload_class_module((str, nuclei.nuclei))
    def add_nuclei(self, label, nuclei_input):
        if label in self.nuclei.keys():
            raise ValueError(f'{label} is already present in the spin system.')
            
        self.nuclei[label] = nuclei_input
        self.coupling_constants.add_nuclei(label)
        self.basis = Basis.Basis(self.nuclei, self.basis.basis_type)
        
        self.name = self.__name_spin_system(self.nuclei)
        self.dir_out_name = self.__make_out_dir(self.name)
        print('New spin sytem set')
        print(f'Outputs will be placed in directory {self.dir_out_name}')
        
        Text.write_spinsystem(self)

            
    def spin_permutation(self, permutations):
        """
        defines spin permutations.

        Parameters
        ----------
        permutations : TYPE: list/array
            DESCRIPTION: contains equivalent nuclei.

        Raises
        ------
        SyntaxError
            DESCRIPTION: raised if permutations contains less than 2 nuclei.
        ValueError
            DESCRIPTION: raised if permutations contains different isotope types
            of if permutations contains a nuclei already present in another group.

        Returns
        -------
        None.

        """
        
        exceptions.type_check(permutations, (list, np.ndarray))
        if len(permutations) < 2:
            raise SyntaxError('At least two spins are required to define a spin permutation.')
        exceptions.duplicate_check(permutations)
        for nuclei_in in permutations:
            exceptions.nuclei_is_defined(list(self.nuclei.keys()), nuclei_in)
            
        isotope_ref = self.nuclei[permutations[0]].isotope
        for nuclei_in in permutations[1:]:
            if self.nuclei[nuclei_in].isotope != isotope_ref:
                raise ValueError('Spin permutation can only be defined for identical isotopes.')
            
        if self.equivalent_nuclei == None:
            self.equivalent_nuclei = dict()
            self.equivalent_nuclei['group1'] = np.asarray(permutations)
            
        else:
            equivalent_nuclei = np.concatenate(list(self.equivalent_nuclei.values()))
            for nuclei_in in permutations:
                if nuclei_in in equivalent_nuclei:
                    raise ValueError(f'Cannot include spin {nuclei_in} in this new group of equivalent nuclei as {nuclei_in} is already present in another group. Use reinitialize_permutation.')
            
            self.equivalent_nuclei[f'group{len(self.EquivalentNuclei)+1}'] = np.asarray(permutations)
            
            
    def reinitialize_permutation(self):
        """
        reinitializes the symmetry.

        Returns
        -------
        None.

        """
        
        self.equivalent_nuclei = None
        
        
    def show_spin_system(self, save = False, figname = None, overwrite = False):
        """
        plots the spin system in a 3d axis frame

        Parameters
        ----------
        save : TYPE: bool, optional
            DESCRIPTION: whether the figure is saved or not. The default is False.
        figname : TYPE: str, optional
            DESCRIPTION: figure name. The default is None.
        overwrite : TYPE: bool, optional
            DESCRIPTION: whether overwriting is allowed or not. The default is False.

        Raises
        ------
        AssertionError
            DESCRIPTION: raised if overwrite is True and figname already exists.

        Returns
        -------
        None.

        """
        
        exceptions.type_check(save, bool)
        
        if save:
            exceptions.type_check(figname, str)
            
            exceptions.folder_create(f'{self.output_directory}/Figures')
            exceptions.folder_create(f'{self.output_directory}/Figures/SpinSystem')
            
            if os.path.isfile(f'{self.output_directory}/Figures/SpinSystem/{figname}'):
                exceptions.type_check(overwrite, bool)
                if overwrite == True:
                    Figures.spin_system_geometry(self.nuclei, save = True, figname = '{self.output_directory}/{figname}')
                else:
                    raise AssertionError('{figname} already exists. Cannot replace.\
                                         Use the option overwrite=True')
        else:
            Figures.spin_system_geometry(self.nuclei, save = False)