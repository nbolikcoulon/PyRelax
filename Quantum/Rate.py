###This modules deals with hamiltonians for spin interactions
#Hamiltonian are given in units of B0
### 

import sys, os

import numpy as np

sys.path.insert(1, '../NMR')

import constants as cst
import SpinOperators
import Interactions





def GetpRange(Type):
    
    if Type == 'DD':
        pRanges = {'2': [0],
                '1': [0, 1],
                '0':[-1, 0, 1],
                '-1': [0, 1],
                '-2': [0]}
    
        return pRanges
        
    elif Type == 'Q':
        pRanges = {'2': [0],
                '1': [0],
                '0':[0],
                '-1': [0],
                '-2': [0]}
        
        return pRanges

    elif 'CSA' in Type:
        pRanges = {'2': [0],
                '1': [0],
                '0':[0],
                '-1': [0],
                '-2': [0]}
        
        return pRanges


    else:
        print('Unknown error occured. Please report to us.')
        
        
        
def GetCorrelationInfo(correlation):
    
    C1, C2 = correlation.split('/')
    
    try:
        T1, N1 = C1.split('-')
    except ValueError:
        T1, N1a, N1b = C1.split('-')
        N1 = [N1a, N1b]
        
    try:
        T2, N2 = C2.split('-')
    except ValueError:
        T2, N2a, N2b = C2.split('-')
        N2 = [N2a, N2b]
        
    return T1, T2, N1, N2



#calculate interaction amplitude
def GetAmplitude(SpinSystem, IntType, Nuclei):
    
    if IntType == 'DD':
        nuclei1, nuclei2 = Nuclei
        
        dd_12 = cst.DipolarConstant(SpinSystem, nuclei1, nuclei2)
        
        return dd_12
    
    elif IntType == 'CSA':
        
        csa = SpinSystem.CSA[Nuclei]['axial']
        
        return csa
    
    elif IntType == 'CSA1':
        
        csa = SpinSystem.CSA[Nuclei]['sig1']
        
        return csa
        
    elif IntType == 'CSA2':
        
        csa = SpinSystem.CSA[Nuclei]['sig2']
        
        return csa
        
    elif IntType == 'Q':
        
        ms = cst.SpinQuantumNumber(SpinSystem, Nuclei)
        
        quad = 168e3 * 2. * np.pi /(ms * (2.*ms - 1))
        
        return quad
    
    
    
#get eigenfrequencies
def _getFrequency(SpinSystem, w):
    
    w_init = w.split('w')
    
    if w_init[0] != '':
        Scalar = float(w_init[0])
        
    else:
        Scalar = 1
        
        
    Nuclei = w_init[-1]
    Isotope = SpinSystem.Nuclei[Nuclei]
    
    return Scalar * cst.GyromagneticRatios[Isotope]
    
    
    
    
def GetFrequency(SpinSystem, w):
    
    if w == '0':
        return 0.
    
    
    #check if the frequency is a sum
    w_sum = w.split('+')
    if len(w_sum) == 2:
        w1, w2 = w_sum
        
        w1_scal = _getFrequency(SpinSystem, w1)
        w2_scal = _getFrequency(SpinSystem, w2)
        
        return w1_scal + w2_scal
    
    #check if the frequency is a difference
    w_diff = w.split('-')
    if len(w_diff) == 2:
        w1, w2 = w_diff
        
        w1_scal = _getFrequency(SpinSystem, w1)
        w2_scal = _getFrequency(SpinSystem, w2)
        
        return w1_scal - w2_scal
    
    #none of the above
    w_scal = _getFrequency(SpinSystem, w)
    
    return w_scal




class Rate:
    
    """"
    compute relaxation rates, as well as associated functions
    standard computing/ploting functions
    use magnetic equivalence to potentially symplify rates
    output formating in python and latex
    
    initializating takes a spin system class and two (not necessarily normalized) spin operators
    all contribution to relaxation are computed. See the Interaction module to obtain interaction-specific relaxation
    """



    ##########################################################
    #                                                        #
    #                      Compute Rate                      #
    #                                                        #
    ##########################################################

    #initialization: computation of the double commutator
    def __init__(self, SpinSystem, Op1, Op2, SecularApproximation = True):
        
        #check dimensions
        Dim = SpinOperators.HilbertDimension(SpinSystem)
        
        if len(Op1) != Dim:
            print("Operator 1 dimension doesn't match Liouville space dimension")
            return None
        
        if len(Op2) != Dim:
            print("Operator 2 dimension doesn't match Liouville space dimension")
            return None


        #initialization
         #Normalization
        Op1N = SpinOperators.NormalizeOperator(SpinSystem, Op1)
        Op2N = SpinOperators.NormalizeOperator(SpinSystem, Op2)
        
        Nuclei = list(SpinSystem.Nuclei.keys())
        
        InteractionsList = dict()
         #dipole-dipole
        for nuclei1 in range(len(Nuclei)-1):
            for nuclei2 in range(nuclei1+1, len(Nuclei)):
                InteractionsList['DD-' + Nuclei[nuclei1] + '-' +  Nuclei[nuclei2]] = Interactions.DDtensors(SpinSystem, Nuclei[nuclei1], Nuclei[nuclei2])
         #CSA
        for nuclei in Nuclei:
            if list(SpinSystem.CSA[nuclei].keys()) == ['axial']:
                InteractionsList['CSA-' + nuclei] = Interactions.CSAtensors(SpinSystem, nuclei)
            elif list(SpinSystem.CSA[nuclei].keys()) == ['sig1', 'sig2']:
                InteractionsList['CSA1-' + nuclei] = Interactions.CSAtensors(SpinSystem, nuclei)
                InteractionsList['CSA2-' + nuclei] = Interactions.CSAtensors(SpinSystem, nuclei)
                
         #quadrupolar
        for nuclei in Nuclei:
            if cst.SpinQuantumNumber(SpinSystem, nuclei) != '1/2':
                InteractionsList['Q-' + nuclei] = Interactions.Qtensors(SpinSystem, nuclei)
        




        #Calculation
        self.DoubleComm = dict()
        
        InteractionsList_types = list(InteractionsList.keys())
        
        for Int_i_count, Int_i in enumerate(InteractionsList_types):
            for Int_j in InteractionsList_types[Int_i_count:]:
                
                Tensors_i, Freq_i = InteractionsList[Int_i]
                Tensors_j, Freq_j = InteractionsList[Int_j]
                
                pRanges_i = GetpRange(Int_i.split('-')[0])
                pRanges_j = GetpRange(Int_j.split('-')[0])
                
                rate = Int_i + '/' + Int_j
                rate_ij = dict()
                
                for q in range(-2, 3):
                    for p1 in pRanges_i[str(q)]:
                        for p2 in pRanges_j[str(q)]:
                            
                            if SecularApproximation:
                                
                                if Freq_i[str(q) + ', ' + str(p1)][1] == Freq_j[str(q) + ', ' + str(p2)][1]:
                                
                                    DoubleCommutator = SpinOperators.Commutator( Tensors_i[str(q) + ', ' + str(p1)], 
                                                                                SpinOperators.Commutator(np.conj(Tensors_j[str(q) + ', ' + str(p2)]).T,
                                                                                                         Op2N)
                                                                                )
                                    
                                    ScalarProduct = SpinOperators.ScalarProduct(Op1N, DoubleCommutator)
                                    Trace = np.trace(ScalarProduct)
                                
                                    if Trace != 0:
                                        if Freq_j[str(q) + ', ' + str(p2)][0] in rate_ij:
                                            try:
                                                rate_ij[Freq_j[str(q) + ', ' + str(p2)][0]]['0, ' + str(q)] += Trace
                                            except KeyError:
                                                rate_ij[Freq_j[str(q) + ', ' + str(p2)][0]]['0, ' + str(-q)] += Trace
                                            
                                        elif str(q) + ', ' + str(-p2) in Freq_j.keys() and Freq_j[str(q) + ', ' + str(-p2)][0] in rate_ij:
                                            try:
                                                rate_ij[Freq_j[str(q) + ', ' + str(-p2)][0]]['0, ' + str(q)] += Trace
                                            except KeyError:
                                                rate_ij[Freq_j[str(q) + ', ' + str(-p2)][0]]['0, ' + str(-q)] += Trace
                                            
                                        else:
                                            rate_ij[Freq_j[str(q) + ', ' + str(p2)][0]] = {'0, ' + str(q): Trace}
                                            
                                            
                                
                                
                            elif not SecularApproximation:
            
                                DoubleCommutator = SpinOperators.Commutator( Tensors_i[str(q) + ', ' + str(p1)], 
                                                                            SpinOperators.Commutator(np.conj(Tensors_j[str(q) + ', ' + str(p2)]).T,
                                                                                                     Op2N)
                                                                            )
                                
                                ScalarProduct = SpinOperators.ScalarProduct(Op1N, DoubleCommutator)
                                Trace = np.trace(ScalarProduct)
            
                
                                if Trace != 0:
                                    
                                    DeltaFreq = Freq_i[str(q) + ', ' + str(p1)][0] + '-' + Freq_j[str(q) + ', ' + str(p2)][0]       #frequency for the oscillating part
                                    key = DeltaFreq + ', ' + str(q)
                                    
                                    try:
                                        rate_ij[Freq_j[str(q) + ', ' + str(p1)][0]][key] += Trace     #case where the frequency of the spectral density and the frequency for the oscillating part have already been encountered
                                        
                                    except:
                                        rate_ij[Freq_j[str(q) + ', ' + str(p1)][0]] = {key: Trace}
            
            
                if rate_ij != {}:
                    self.DoubleComm[rate] = rate_ij
                    
                    
                    
        
                
                
                
            
                    
                    
                    
                    
    #compute relaxation rate
    def EvalRate(self, SpinSystem, Model, B0, ParametersDynamics, TimeDpdtInteractionsAmp = False, SecularApproximation = True, time = 0):
        
        args = ParametersDynamics
        
        CorrelationsList = list(self.DoubleComm.keys())
        
        rate = 0.
        
        #loop on all interactions -> get the product of the interactions amplitudes and then eval the spectral density at the frequencies of the double commutators
        
        for Corr in CorrelationsList:
            
            Type1, Type2, Nuclei1, Nuclei2 = GetCorrelationInfo(Corr)
            
            amplitude1 = 1./2. * GetAmplitude(SpinSystem, Type1, Nuclei1)   #Factor 1/2 comes from the definition of the correlation function/spectral density function
            amplitude2 = 1./2. * GetAmplitude(SpinSystem, Type2, Nuclei2)
            
            if Type1 == 'DD':
                amplitude1 = {'-2': amplitude1 * np.sqrt(6.),
                              '-1': amplitude1 * np.sqrt(6.),
                              '0': amplitude1 * np.sqrt(6.),
                              '1': amplitude1 * np.sqrt(6.),
                              '2': amplitude1 * np.sqrt(6.)}
                
            elif 'CSA' in Type1:
                amplitude1 = {'-2': amplitude1 * np.sqrt(2./3.) * B0 * cst.GyromagneticRatios[SpinSystem.Nuclei[Nuclei1]],
                              '-1': amplitude1 * np.sqrt(2./3.) * B0 * cst.GyromagneticRatios[SpinSystem.Nuclei[Nuclei1]],
                              '0': amplitude1 * np.sqrt(2./3.) * B0 * cst.GyromagneticRatios[SpinSystem.Nuclei[Nuclei1]],
                              '1': amplitude1 * np.sqrt(2./3.) * B0 * cst.GyromagneticRatios[SpinSystem.Nuclei[Nuclei1]],
                              '2': amplitude1 * np.sqrt(2./3.) * B0 * cst.GyromagneticRatios[SpinSystem.Nuclei[Nuclei1]]}
                
            elif Type1 == 'Q':
                amplitude1 = {'-2': amplitude1 * SpinSystem.EFG_asymmetry[Nuclei1],
                              '-1': 0.,
                              '0': amplitude1 * np.sqrt(6.),
                              '1': 0.,
                              '2': amplitude1 * SpinSystem.EFG_asymmetry[Nuclei1]}
                
                
            if Type2 == 'DD':
                amplitude2 = {'-2': amplitude2 * np.sqrt(6.),
                              '-1': amplitude2 * np.sqrt(6.),
                              '0': amplitude2 * np.sqrt(6.),
                              '1': amplitude2 * np.sqrt(6.),
                              '2': amplitude2 * np.sqrt(6.)}
                
            elif 'CSA' in Type2:
                amplitude2 = {'-2': amplitude2 * np.sqrt(2./3.) * B0 * cst.GyromagneticRatios[SpinSystem.Nuclei[Nuclei2]],
                              '-1': amplitude2 * np.sqrt(2./3.) * B0 * cst.GyromagneticRatios[SpinSystem.Nuclei[Nuclei2]],
                              '0': amplitude2 * np.sqrt(2./3.) * B0 * cst.GyromagneticRatios[SpinSystem.Nuclei[Nuclei2]],
                              '1': amplitude2 * np.sqrt(2./3.) * B0 * cst.GyromagneticRatios[SpinSystem.Nuclei[Nuclei2]],
                              '2': amplitude2 * np.sqrt(2./3.) * B0 * cst.GyromagneticRatios[SpinSystem.Nuclei[Nuclei2]]}
                
            elif Type2 == 'Q':
                amplitude2 = {'-2': amplitude2 * SpinSystem.EFG_asymmetry[Nuclei2],
                              '-1': 0.,
                              '0': amplitude2 * np.sqrt(6.),
                              '1': 0.,
                              '2': amplitude2 * SpinSystem.EFG_asymmetry[Nuclei2]}
                
                
                
            
            for Freq in self.DoubleComm[Corr].keys():
                
                w_density = B0 * GetFrequency(SpinSystem, Freq)
                w_secular, q = self.DoubleComm[Corr][Freq][0].split(',')
                w_secular = B0 * GetFrequency(SpinSystem, w_secular)
                
                if not TimeDpdtInteractionsAmp:
                    if SecularApproximation:
                        rate += 1./2. * amplitude1[q] * amplitude2[q] * self.DoubleComm[Corr][Freq][1] * Model(w_density, args, Type1, Type2)
                        
                    if not SecularApproximation:
                        rate += 1./2. * amplitude1[q] * amplitude2[q] * np.exp(1j * w_secular * time) * self.DoubleComm[Corr][Freq][1] * Model(w_density, args, Type1, Type2)
                    
                if TimeDpdtInteractionsAmp:
                    if SecularApproximation:
                        rate += 1./2. * self.DoubleComm[Corr][Freq][1] * Model(w_density, args, Type1, Type2, amplitude1[q], amplitude2[q])
                        
                    if not SecularApproximation:
                        rate += 1./2. * np.exp(1j * w_secular * time) * self.DoubleComm[Corr][Freq][1] * Model(w_density, args, Type1, Type2, amplitude1[q], amplitude2[q])
                        
                        
                
                
            
        
        
        
        
        
        return rate
            
            
            
            
            
            
            
            
            
            
            
            


