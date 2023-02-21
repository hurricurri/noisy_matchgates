import numpy as np 
from noise_channels import global_depolarizing_channel, amplitude_damping, X_bit_flip, global_x_rotation
import random

def list_to_array(lst):
    return np.array(lst)

def quantum_process(rho, U, n, p, angle, noisechannel):
    '''Takes as input an initial state, a unitary to be applied, the size of the system 'n' and
    a probability associated with a specific noise channel, the noise channel and an angle associated to the X rotation (in radians), 
    and measures in the computational basis, outputting a bitstring b of 0s and 1s'''

    rho = U @ (rho @ U.conj().T)
    

    ##apply noise channel based on input
    if noisechannel == 'depolarizing':

        rho = global_depolarizing_channel(n, rho, p) 

    elif noisechannel == 'amplitude damping':

        rho = amplitude_damping(n, rho, p)
    
    elif noisechannel == 'bit flip':

        rho = X_bit_flip(n, rho, p)
    
    elif noisechannel == 'X rotation':

        rho = global_x_rotation(n, rho, p, angle)
    
    else: 
        rho = rho


    ##measure in the computational basis 
    diag = np.diag(rho) 

    probabilities = np.abs(diag)**2
    
    probabilities = probabilities/np.sum(probabilities)

    outputs = [random.choices([0, 1], weights = [1 - p, p])[0] for p in probabilities]

    barray = list_to_array(outputs)
    
    return barray




##Example of how to run:
#n=4
#purestate = pure_initial(n)
#covar = generate_covariance(n)
#c, a = fermionic_operators(n)
#gaussianstate = gaussian_density_matrix(covar, c, a)
#U = random_FGUclifford(n,ret_unitary=True)[1]
#p=0.1
#
#b = quantum_process(gaussianstate, U, n, p, False, 'amplitude damping')
#
#b = quantum_process(purestate, U, n, p, False, 'amplitude damping')
#
#print('b array for pure state:', b)
