import numpy as np 
import math as mt 
import random
from sample_matchgates import X, I 



##Global depolarizing noise channel
def global_depolarizing_channel(n, rho, p):

    max_mix_state= np.identity(2**n, dtype = 'complex128') / 2**n

    dep_state = (1-p)*rho + p*max_mix_state
    
    return dep_state


##Global amplitude damping channel
def amplitude_damping(n, density_matrix, gamma):

    I = np.identity(2, dtype='complex128')

    K0 = np.array([[1, 0], [0, np.sqrt(1 - gamma)]], dtype='complex128')

    K1 = np.array([[0, np.sqrt(gamma)], [0, 0]], dtype='complex128')

    for j in range(n-1):
        K0 = np.kron(K0, I)
        K1 = np.kron(K1, I) 

    ampstate = K0 @ density_matrix @ K0.conj().T + K1 @ density_matrix @ K1.conj().T

    return ampstate


##X bit-flip on all qubits with probability p
def X_bit_flip(n, density_matrix, p):

    bit_flip = (1-p) * I + p * X

    for j in range(n-1):

        bit_flip = np.kron(bit_flip, I)

    flipped_state = bit_flip @ density_matrix @ bit_flip.conj().T

    return flipped_state

##X-rotation on all qubits with angle theta and probability p
def global_x_rotation(n, density_matrix, p, theta):
    
    prob = np.random.randint(1)

    Identity = np.identity(2**n, dtype='complex128')

    rotation = np.array([[np.cos(theta/2), -1j*np.sin(theta/2)], [-1j*np.sin(theta/2), np.cos(theta/2)]], dtype='complex128')
    
    if prob > p: 
        for i in range(n-1):
            rotation = np.kron(rotation, I)
    else: 
            rotation = Identity
    
    rotated_state = rotation @ density_matrix @ rotation.conj().T

    return rotated_state

