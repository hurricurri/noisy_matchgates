import numpy as np

# required gate definitions:
# single-qubit Z rotation

def Z(theta):
    return np.array([[np.exp(-1j * theta / 2), 0],
                     [0, np.exp(1j * theta / 2)]], dtype='complex128')
    
# two-qubit XX rotation
def XX(theta):
    return np.array([[np.cos(theta/2), 0, 0, -1j * np.sin(theta/2)],
                     [0, np.cos(theta/2), -1j * np.sin(theta/2), 0],
                     [0, -1j * np.sin(theta/2), np.cos(theta/2), 0],
                     [-1j * np.sin(theta/2), 0, 0, np.cos(theta/2)]], dtype='complex128')
    
# X gate
X = np.array([[0, 1],
              [1, 0]], dtype='complex128')
              
# Y gate
Y = np.array([[0, -1j],
              [1j, 0]], dtype='complex128')

# single-qubit identity
I = np.identity(2, dtype='complex128')
    


def random_transp(n):
    """returns a random permutation of the numpy array [1,...,n], and also a decomposition into nearest-neighbour transpositions (list of 2-tuples)"""
    
    transp = np.arange(1,n+1, dtype='int')
    nn_transp = []
    
    # Fisher–Yates shuffle
    for i in range(1,n):
        j = np.random.randint(i,n+1)
        
        mem = transp[i-1]
        
        transp[i-1] = transp[j-1]
        transp[j-1] = mem
    
        for k in range(j-i):
            nn_transp.append((i+k,i+k+1))
        
        for k in range(j-i-1):
            nn_transp.append((j-2-k,j-1-k))
    
    transp = (np.argsort(transp)+1)
    
    return transp, nn_transp


def random_FGUclifford(n, ret_unitary = False):
    """outputs uniformly random 2n*2n signed permutation matrix i.e. matchgate Clifford
       if ret_unitary is set to True, the n-qubit unitary is returned as well"""
    
    transp, nn_transp = random_transp(2*n)
    
    Q = np.identity(2*n, dtype='float32')[:,transp-1]
    if ret_unitary == True:
    
        unitary = np.identity(2**n, dtype='complex128')
        
        for i,j in enumerate(nn_transp):

            if j[0] % 2 == 1: # Z-rotation by pi/2 on qubit j[0]//2, identity on all others (up to sign)
    
                if j[0]//2 == 0:
                    gate = Z(np.pi/2)
                    
                    for qubit in range(1,n):
                        gate = np.kron(gate, I)
                       
                else:
                    gate = I
                    
                    for qubit in range(1, n):
                            if qubit == j[0]//2:
                                gate = np.kron(gate, Z(np.pi/2))
                                
                            else:
                                gate = np.kron(gate, I)
                
                unitary = gate @ unitary
                
                # correct sign (for post-processing purposes)
                # apply X to qubit j[1]//2-1, Z to all subsequent qubits
                    
                if j[1]//2 - 1 == 0:
                    gate = X
                        
                    for qubit in range(1,n):
                        gate = np.kron(gate, Z(np.pi))
                    
                else:
                    gate = I
                        
                    for qubit in range(1,j[1]//2-1):
                        gate = np.kron(gate, I)
                        
                    gate = np.kron(gate, X)
                        
                    for qubit in range(j[1]//2,n):
                        gate = np.kron(gate, Z(np.pi))
                        
                unitary = gate @ unitary
                
                
                
            else: # XX-rotation by pi/2 on qubits j[0]//2-1 and j[0]//2, identity on all others (up to sign, but we randomize over all of them afterwards)

                if j[0]//2 - 1 == 0:
                    gate = XX(np.pi/2)
                    
                    for qubit in range(2,n):
                        gate = np.kron(gate,I)
                       
                else:
                    gate = I
                    
                    for qubit in range(1, n):
                            if qubit == j[0]//2 - 1:
                                gate = np.kron(gate, XX(np.pi/2))
                                
                            elif qubit == j[0]//2:
                                continue
                                
                            else:
                                gate = np.kron(gate, I)
                
                unitary = gate @ unitary
                
                # correct sign (for post-processing purposes)
                # apply Y to qubit j[1]//2, Z to all subsequent qubits
                    
                if j[1]//2 == 0:
                    gate = Y
                        
                    for qubit in range(1,n):
                        gate = np.kron(gate, Z(np.pi))
                    
                else:
                    gate = I
                        
                    for qubit in range(1,j[1]//2):
                        gate = np.kron(gate, I)
                        
                    gate = np.kron(gate, Y)
                        
                    for qubit in range(j[1]//2+1,n):
                        gate = np.kron(gate, Z(np.pi))
                  
                    
                    
                unitary = gate @ unitary
                     
    # signed permutations
    for i in range(1,2*n+1):
        if np.random.randint(2) == 1: # flip sign
            
            Q[:,i-1] *= (-1)
            
            if ret_unitary == True:
            
                if i % 2 == 1: # apply Y to qubit i//2, Z to all subsequent qubits
                    
                    if i//2 == 0:
                        gate = Y
                        
                        for qubit in range(1,n):
                            gate = np.kron(gate, Z(np.pi))
                    
                    else:
                        gate = I
                        
                        for qubit in range(1,i//2):
                                gate = np.kron(gate, I)
                        
                        gate = np.kron(gate, Y)
                        
                        for qubit in range(i//2+1,n):
                            gate = np.kron(gate, Z(np.pi))
                  
                    

                else: # apply X to qubit i//2-1, Z to all subsequent qubits
                    
                    if i//2 - 1 == 0:
                        gate = X
                        
                        for qubit in range(1,n):
                            gate = np.kron(gate, Z(np.pi))
                    
                    else:
                        gate = I
                        
                        for qubit in range(1,i//2-1):
                            gate = np.kron(gate, I)
                        
                        gate = np.kron(gate, X)
                        
                        for qubit in range(i//2,n):
                            gate = np.kron(gate, Z(np.pi))
                        
                unitary = gate @ unitary
    
    if ret_unitary == True:
        return Q, np.conj(unitary).T   
        
    return Q

