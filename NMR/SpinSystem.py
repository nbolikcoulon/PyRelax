##This module is a class for setting the spin system. This includes:
# setting the spin nuclei
# setting the geometry (optional)
# setting the CSA properties (optional)
### 
import sys, os

import numpy as np
from fractions import Fraction

import constants as cst 

sys.path.insert(1, '../Outputs')
sys.path.insert(2, '../Quantum')

import Figures
import Text
import SpinOperators
import Interactions
import Rate


####################################################
#                                                  #
#                      Naming                      #
#                                                  #
####################################################

def NameSpinSystem(SpinSystem):
    
    Counts = {'1H': 0, '2H': 0, '15N': 0, '13C': 0, '19F': 0, '31P': 0}
    
    for Isotopes in SpinSystem.values():
        Counts[Isotopes] += 1
        
    Name = ''
    for Isotope, count in Counts.items():
        if count == 0:
            pass
        elif count == 1:
            Name += Isotope + '_'
        else:
            Name += Isotope + str(count) + '_'
            
    Name = Name[:-1]
    
    return Name



def MakeOutDir(name):
    OutDir = os.getcwd() + '/' + name
    
    if not os.path.isdir(OutDir):
        os.mkdir(OutDir)
        
        return OutDir + '/'
    
    else:
        n = 1
        OutDir_new = OutDir + '_' + str(n)
        while os.path.isdir(OutDir_new):
            n += 1
            OutDir_new = name + '_' + str(n)
            
        os.mkdir(OutDir_new)
        
        return OutDir_new + '/'





####################################################
#                                                  #
#                      Checks                      #
#                                                  #
####################################################

#check if input label is in the spin system
def CheckDefinedLabel(SpinSystem, Label):
    if Label not in SpinSystem.keys():
        print('Please enter a label as defined in the spin system: ', str(list(SpinSystem.keys()))[1:-1])
        return False
    else:
        return True


#check if output folder exists
def CheckFolder(outputFolder, Folder):
    if not os.path.isdir(outputFolder + Folder):
        os.mkdir(outputFolder + Folder)
        
        









class SetSpinSystem:
    
    ################################################################
    #                                                              #
    #                      Define spin system                      #
    #                                                              #
    ################################################################
    
    #initialyse the spin system
    def __init__(self, SpinSys):    #input must be a dictionanary of label: isotope, all strings
    
        self.Nuclei = dict()
    
        try:

            for Labels, Isotopes in SpinSys.items():
                try:
                    cst.SpinQuantumNumber(self, Labels)
                    self.Nuclei[str(Labels)] = Isotopes
                    
                    if cst.SpinQuantumNumber(self, Labels) != '1/2':
                        self.EFG_asymmetry[str(Labels)] = 0.        #asymmetry of the EFG quadrupolar tensor, defined as abs{(Vyy - Vxx)/Vzz}
                        
                    
                except KeyError:
                    print(f'The nuclei {Labels} cannot be found in our list of possible isotopes. Please make sure it exists or contact us to add it to our list.')
                    print('Spin system not defined.')
                    self.Nuclei.clear()
                    break
                
            self.CSA = dict()
            self.Coordinates = dict()
            for count, Label in enumerate(SpinSys.keys()):
                self.CSA[Label] = {'isotropic': 0.0}
                self.Coordinates[Label] = []
            self.BasisType = None
            self.Basis = None    
            self.CouplingConstants = np.zeros( (len(SpinSys), len(SpinSys)) )
            self.EquivalentNuclei = None
                
            print('Spin sytem set')
            
            self.Name = NameSpinSystem(self.Nuclei)
            self.DirOutName = MakeOutDir(self.Name)
            print(f'Outputs will be placed in directory {self.DirOutName[:-1]}')
            
            Text.write_spinsystem(self, self.DirOutName)
        
            print()
            
            
        except AttributeError:
            print("The spin system must be set using a dictionnary: {'Label': 'Isotope'}.")
            
            
    #remove one isotope from spin system        
    def remove_Nuclei(self, Label):
        try:
            
            Position = np.where(np.asarray(list(self.Nuclei.keys)) == Label)[0][0]
            
            del self.Nuclei[Label]
            del self.CSA[Label]
            del self.Coordinates[Label]
            
            CouplingSaved = self.CouplingConstants
            self.CouplingConstants = np.zeros( (len(self.Nuclei), len(self.Nuclei)) )
            self.CouplingConstants[:Position, :Position] = CouplingSaved[:Position, :Position]
            self.CouplingConstants[Position:, Position:] = CouplingSaved[Position+1:, Position+1:]

            
            self.Name = NameSpinSystem(self.Nuclei)
            self.DirOutName = MakeOutDir(self.Name)
            print('New spin sytem set')
            print(f'Outputs will be placed in directory {self.DirOutName}')
            
            Text.write_spinsystem(self, self.DirOutName)

        except KeyError:
            print(f'{Label} is not present in the spin system')
            
    #add one nuclei in the spin system
    def add_nuclei(self, Label, Isotope):
        try:
            cst.SpinQuantumNumber(self, Label)
            self.Nuclei[str(Label)] = Isotope
            self.CSA[Label] = {'isotropic': 0.0}
            self.Coordinates[Label] = []
            
            CouplingSaved = self.CouplingConstants
            self.CouplingConstants = np.zeros( (len(self.Nuclei), len(self.Nuclei)) )
            self.CouplingConstants[:-1, :-1] = CouplingSaved

            
            
            print('New spin sytem set')
            
            self.Name = NameSpinSystem(self.Nuclei)
            self.DirOutName = MakeOutDir(self.Name)
            print(f'Outputs will be placed in directory {self.DirOutName[:-1]}')
            
            Text.write_spinsystem(self, self.DirOutName)
        
            print()
            
        except KeyError:
            print(f'The nuclei {Label} cannot be found in our list of possible isotopes. Please make sure it exist or contact us to add it to our list.')
            print('Nuclei not added')
            
            
            
            
            
    #############################################################
    #                                                           #
    #                      Scalar coupling                      #
    #                                                           #
    #############################################################
    
    def DefineSalarCoupling(self, Nuclei1, Nuclei2, J12):
        if CheckDefinedLabel(self.Nuclei, Nuclei1):
            if CheckDefinedLabel(self.Nuclei, Nuclei2):

                Position1 = np.where(np.asarray(list(self.Nuclei.keys())) == Nuclei1)[0][0]
                Position2 = np.where(np.asarray(list(self.Nuclei.keys())) == Nuclei2)[0][0]
                
                try:
                    self.CouplingConstants[Position1, Position2] = float(J12)
                    self.CouplingConstants[Position2, Position1] = float(J12)
                    
                    Text.write_spinsystem(self, self.DirOutName)
                    
                except ValueError:
                    print('The scalar coupling constant must be a scalar')


            
            
            
            
    ######################################################
    #                                                    #
    #                      Geometry                      #
    #                                                    #
    ######################################################
            
    #add coordinates for one nuclei. Can be Cartesian [x, y, z] or Polar [theta, phi, r] with angles in rad or deg
    def add_Coordinate(self, Label, Vector, system = 'Cartesian', angles = 'rad'):
        if CheckDefinedLabel(self.Nuclei, Label):
            try:
                if system == 'Cartesian':
                    x, y, z = Vector
                    self.Coordinates[Label] = np.asarray([float(x), float(y), float(z)])
                    
                    print(f'Coordinates of Nuclei {Label} set to x = {self.Coordinates[Label][0]}, y = {self.Coordinates[Label][1]}, z = {self.Coordinates[Label][2]}')
                    print()
                    
                    Text.write_spinsystem(self, self.DirOutName)
 
                elif system == 'Polar':
                    
                    theta, phi, r = Vector
                    
                    if angles == 'rad':
                        theta, phi, r = float(theta), float(phi), float(r)
                        x, y, z = r * np.sin(theta) * np.cos(phi), r * np.sin(theta) * np.sin(phi), r * np.cos(theta)
                        self.Coordinates[Label] = np.asarray([x, y, z])
                        
                        print(f'Coordinates of Nuclei {Label} set to x = {self.Coordinates[Label][0]}, y = {self.Coordinates[Label][1]}, z = {self.Coordinates[Label][2]}')
                        print()
                        
                        Text.write_spinsystem(self, self.DirOutName)
                        
                    elif angles == 'deg':
                        theta, phi, r = np.deg2rad(float(theta)), np.deg2rad(float(phi)), np.deg2rad(float(r))
                        x, y, z = r * np.sin(theta) * np.cos(phi), r * np.sin(theta) * np.sin(phi), r * np.cos(theta)
                        self.Coordinates[Label] = np.asarray([x, y, z])
                        
                        print(f'Coordinates of Nuclei {Label} set to x = {self.Coordinates[Label][0]}, y = {self.Coordinates[Label][1]}, z = {self.Coordinates[Label][2]}')
                        print()
                        
                        Text.write_spinsystem(self, self.DirOutName)
                            
                    else:
                        print('The angle unit must be either rad or deg. Default is rad.')
    
                else:
                    print('Coordinate system frame must be either Cartesian or Polar. Default is Cartesian.')
                    
            except ValueError:
                print('Coordinates must be for a 3D axis system.')
            

    #change the coordinate for one nuclei
    def change_Coordinate(self, Label, Vector, system = 'Cartesian', angles = 'rad'):
        self.add_Coordinate(Label, Vector, system, angles)

        
    #remove the coordinate for one nuclei
    def remove_Coordinate(self, Label):
        if CheckDefinedLabel(self.Nuclei, Label):
            self.Coordinates[Label] = []
            print(f'Coordinates of Nuclei {Label} have been removed')
            print()
            
            Text.write_spinsystem(self, self.DirOutName)

            
    #visualize the spin system
    def ShowSpinSystem(self, save = False):
        if save:
            CheckFolder(self.DirOutName, 'Figure')
            CheckFolder(self.DirOutName, 'Figure/Spinsystem')
            
            overwrite = False
            if os.path.isfile(self.DirOutName + 'Figure/Spinsystem/Geometry.pdf'):
                
                choice = str(input('Figure Geometry.pdf already exits. Do you want to overwrite it?  '))
                if choice == 'yes':
                    overwrite = True
                    
            Figures.SpinSystemGeometry(self, self.Coordinates, save = True, outputFolder = self.DirOutName + '/Figure/Spinsystem/', overwrite = overwrite)
        else:
            Figures.SpinSystemGeometry(self, self.Coordinates, save = False)
                
            
            
            
            
    ########################################################
    #                                                      #
    #                      Symmetrize                      #
    #                                                      #
    ########################################################
            
    #set equivalent in the spin system. Permunation is an array containing labels as defined in the spin system
    def SpinPermutation(self, Permutation = None):
        
        #check all spins are defined in the spin system
        for nuclei in Permutation:
            if not CheckDefinedLabel(self.Nuclei, nuclei):
                        
                return False
            
        #check Permutation contains at least two items
        if len(Permutation) < 2:
            print('At least two spins are required to use this functions.')
            
            return False
        
        try:
            #initialize if necessary
            if self.EquivalentNuclei == None:
                    
                self.EquivalentNuclei = dict()
                
                self.EquivalentNuclei['group1'] = Permutation
                
            else:
                
                #check Permutation does not contain nuclei already present in other groups
                EquivalentNuclei = list(self.EquivalentNuclei.values())
                for spin in Permutation:
                    if spin in EquivalentNuclei:
                        print(f'Cannot include spin {spin} in this new group of equivalent nuclei as {spin} is already present in another group.')
                        
                        return False
                
                
                Ngroups = len(list(self.EquivalentNuclei.keys()))
                
                self.EquivalentNuclei[f'group{Ngroups+1}'] = Permutation
            
            
        except:
            print('Unknown error occured. Please report to us.')
        
        
        
        
        
        

    #################################################
    #                                               #
    #                      CSA                      #
    #                                               #
    #################################################        

    #change the CSA for one residue
    def change_CSA(self, Label, symmetry, amplitude):
        if CheckDefinedLabel(self.Nuclei, Label):
            if symmetry not in ['isotropic', 'axial', 'asymmetric']:
                print("The symmetry must be set to 'isotropic', 'axial' or 'asymmetric'. CSA not changed.")
                
            else:
                if symmetry == 'isotropic':
                    try:
                        if isinstance(amplitude, int) or isinstance(amplitude, float):
                            self.CSA[Label] = {'isotropic': float(amplitude)}
                            
                        elif isinstance(amplitude, list) and np.shape(np.asarray(amplitude)) == (1,):
                            self.CSA[Label] = {'isotropic': float(amplitude[0])}
                            
                        elif isinstance(amplitude, np.ndarray) and np.shape(amplitude) == (1,):
                            self.CSA[Label] = {'isotropic': float(amplitude[0])}
                            
                        else:
                            print('CSA value must be set to a floating number. CSA not changed.')
                            
                            
                        print(f'CSA for {Label} set to isotropic with value ', self.CSA[Label]['isotropic'])
                        print()
                        Text.write_spinsystem(self, self.DirOutName)
                            
                    except:
                        print('An unknown error occured. Please report to us.')
                        
                        
                elif symmetry == 'axial':
                    try:
                        if isinstance(amplitude, int) or isinstance(amplitude, float):
                            self.CSA[Label] = {'axial': float(amplitude)}
                            
                        elif isinstance(amplitude, list) and np.shape(np.asarray(amplitude)) == (1,):
                            self.CSA[Label] = {'axial': float(amplitude[0])}
                            
                        elif isinstance(amplitude, np.ndarray) and np.shape(amplitude) == (1,):
                            self.CSA[Label] = {'axial': float(amplitude[0])}
                            
                        else:
                            print('CSA value must be set to a floating number. CSA not changed.')
                            
                        print(f'CSA for {Label} set to axial with value ', self.CSA[Label]['axial'])
                        print()
                        Text.write_spinsystem(self, self.DirOutName)
                            
                    except:
                        print('An unknown error occured. Please report to us.')
                        
                        
                elif symmetry == 'asymmetric':
                    try:
                        sigma1, sigma2 = amplitude
                        self.CSA[Label] = {'sig1': float(sigma1),
                                           'sig2': float(sigma2)}
                        
                        print(f'CSA for {Label} set to asymmetric with values ', self.CSA[Label]['sig1'], ' and ', self.CSA[Label]['sig2'])
                        print()
                        Text.write_spinsystem(self, self.DirOutName)
                            
                    except ValueError:
                        print('Two CSA amplitudes must be provided for asymmetric CSA in a form of an array. CSA not changed.')
                            

    #add CSA for one nuclei
    def add_CSA(self, Label, symmetry, amplitude):
        self.change_CSA(Label, symmetry, amplitude)
        
        
    #remove CSA for one nuclei. It will be set to isotropic with value 0.0
    def remove_CSA(self, Label):
        if CheckDefinedLabel(self.Nuclei, Label):
            self.CSA[Label] = {'isotropic': 0.0}
            print(f'CSA for {Label} has been removed.')
            
            Text.write_spinsystem(self, self.DirOutName)
            
            
            
    ###########################################################
    #                                                         #
    #                      Asymmetry EFG                      #
    #                                                         #
    ###########################################################
    
    def SetEFG(self, Label, Value):
        
        if CheckDefinedLabel(self.Nuclei, Label):
            
            if cst.SpinQuantumNumber(self, Label) == '1/2':
                print('Asymmetry of EFG tensor can only be set for non spin-1/2 nuclei')
                print()
                
            else:
                try:
                    self.EFG_asymmetry[Label] = float(Value)
                    
                except ValueError:
                    print('A valid value for the asymmetry of the EFG tensor has to be provided')
                    print()
            
            
            
            
    ###################################################
    #                                                 #
    #                      Basis                      #
    #                                                 #
    ###################################################
    
    #print the Operator Basis
    def PrintBasis(self):
        
        if self.Basis == None:
            print('Basis not defined yet. Run SetBasis first')
        
        else:
            print(f'Basis set to {self.BasisType}.')
            print(list(self.Basis.keys()))
            print()
            
    #Define a basis
    def SetBasis(self, OperatorType):   #OperatorType is a string
        
        try:
            self.Basis = eval('SpinOperators.' + OperatorType)(self)
            self.BasisType = OperatorType
            
            self.PrintBasis()
            
            Text.write_spinsystem(self, self.DirOutName)
            
        except:
            print('Please use one of the following options to define the basis:')
            print(' - CartesianOperatorBasis')
            print(' - ShiftOperatorBasis')
            print(' - ZeemanBasis')
            print(' - SingletTripletBasis (for two spin-1/2 spin systems)')
            print()
            
            self.BasisType = None
            self.Basis = None
            
            
            
    def SetCustomBasis(self, OperatorSet_in):
        
        OperatorSet = SpinOperators.CheckBasis(self, OperatorSet_in)
        
        if OperatorSet != False:
        
            print('Setting the basis...')
            self.BasisType = 'Custom'
            self.Basis = OperatorSet
            print(' Setting the basis: Done')
                
            Text.write_spinsystem(self, self.DirOutName)
            
        
        
        
            
        
            
            
########################################
########################################
#below are examples of how to use the current version        
        
            
S = dict()
S['N'] = '15N'
S['HA'] = '1H'
S['HB'] = '1H'
A = SetSpinSystem(S)


#define operators
op1 = SpinOperators.SpinOp(A, 'z', 'N', Normalize = False)
op2 = SpinOperators.SpinOp(A, 'z', 'HA', Normalize = False) @ SpinOperators.SpinOp(A, '+', 'N', Normalize = False)
op3 = SpinOperators.SpinOp(A, 'z', 'HA', Normalize = False) @ SpinOperators.SpinOp(A, '-', 'N', Normalize = False)

#scalar product
Scal = SpinOperators.ScalarProduct(op1, op1)

#normalization
op = SpinOperators.NormalizeOperator(A, op3 @ op4)

#commutation
a = SpinOperators.Commutator(SpinOperators.SpinOp(A, '+', 'HA', Normalize = True), SpinOperators.SpinOp(A, '-', 'HA', Normalize = True))

#relaxation rate
R = Rate.Rate(A, op1, op1)

#double commutator part of the rate
print(R.DoubleComm)

#add characteristics to the spin system
A.add_Coordinate('HA', [0, 1, 0], 'Cartesian')
A.add_Coordinate('HB', [0, -1, 0], 'Cartesian')
A.add_Coordinate('N', [-0.7, 0, 0.7], 'Polar')
A.add_CSA('HA', 'isotropic', [9.0])
A.add_CSA('HB', 'axial', 9.0)
A.add_CSA('N', 'axial', 9.0)

A.DefineSalarCoupling('HA', 'N', 1.0)

#customized basis
opE = SpinOperators.SpinOp(A, 'E', 'HA', Normalize = False) 
opp = SpinOperators.SpinOp(A, '+', 'HA', Normalize = False) 
opm = SpinOperators.SpinOp(A, '-', 'HA', Normalize = False) 
opz = SpinOperators.SpinOp(A, 'z', 'HA', Normalize = False) 
CustBasis = {'E': opE, '+': opp, '-': opm, 'z': opz}

A.SetCustomBasis(CustBasis)



#decompose an operator one a basis
op = SpinOperators.SpinOp(A, 'x', 'HA', Normalize = False) 
Decomp = SpinOperators.BasisProjection(A.Basis, op)


#set a new basis
A.SetBasis('ShiftOperatorBasis')

Ket = SpinOperators.Ket(A, 'aab')

#coherence order of spin operator
print(SpinOperators.CoherenceOrder(A, op))

