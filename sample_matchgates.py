import numpy as np
import numpy.linalg

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
    

def sample_sphere(d):
    """Returns the first component of a random point on a d-dimensional sphere"""

    u = np.random.normal(0, 1 ,d)
    
    # normalize and take first component:
    x = u[0] / np.linalg.norm(u)
    
    return x


def random_matchgate(n):
    """Returns a random generalized matchgate unitary on n qubits"""

    # actual construction of the matchgate unitary:
    unitary = np.identity(2**n, dtype='complex128')
    
    for i in range(2*n - 1, 0, -1):
        for j in range(i, 2*n):
            
            angle = np.arccos(sample_sphere(2*n - (j+1) + 2))
            
            if j % 2 == 0:  # even j
                # XX-rotation by given angle on qubits j//2 - 1 and j//2, identity on all others:
                if j//2 - 1 == 0:
                    gate = XX(angle)
                    
                    for qubit in range(2, n):
                        gate = np.kron(gate, I)
                        
                else:
                    gate = I
                    
                    for qubit in range(1, n):
                        if qubit == j//2 - 1:
                            gate = np.kron(gate, XX(angle))
                            
                        elif qubit == j//2:
                            continue
                            
                        else:
                            gate = np.kron(gate, I)
 
            else:  # odd j
                # Z-rotation by given angle on qubit j//2, identity on all others:
                if j//2 == 0:
                    gate = Z(angle)
                    
                    for qubit in range(1, n):
                        gate = np.kron(gate, I)
                        
                else:
                    gate = I
                    for qubit in range(1, n):
                        if qubit == j//2:
                            gate = np.kron(gate, Z(angle))
                        else:
                            gate = np.kron(gate, I)
                        
            unitary = gate @ unitary
    
    
    if np.random.randint(2):
        # X gate on the last qubit, identity on all others:
        if n == 1:
            gate = X
        else:
            gate = I
            for qubit in range(1, n-1):
                gate = np.kron(gate, I)
            gate = np.kron(gate, X)
        
        unitary = gate @ unitary
    
    return unitary




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
    
    return transp, nn_transp


def random_FGUclifford(n, ret_unitary = False):
    """outputs uniformly random 2n*2n signed permutation matrix i.e. matchgate Clifford
       if ret_unitary is set to True, the n-qubit unitary is returned as well"""
    
    transp, nn_transp = random_transp(2*n)
    Q = np.identity(2*n, dtype='float32')[:,transp-1]
    if ret_unitary == True:
    
    	unitary = np.identity(2**n, dtype='complex128')
    	
    	for i,j in enumerate(nn_transp):

		    if j[0] % 2 == 1: # Z-rotation by pi/2 on qubit j[0]//2, identity on all others
		    
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
		        
		    else: # XX-rotation by pi/2 on qubits j[0]//2-1 and j[0]//2, identity on all others

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
                     
    # signed permutations
    for i in range(1,2*n+1):
        if np.random.randint(2) == 1: # flip sign
            
            Q[i-1,:] *= (-1)
            
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
        return Q, unitary   
        
    return Q
