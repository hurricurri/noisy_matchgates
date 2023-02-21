import numpy as np
import itertools
import scipy.special
import matplotlib.pyplot as plt

from sample_matchgates import random_FGUclifford, I, X, Y
from quantum_process import quantum_process
from fermionic_states import pure_initial, gaussian_density_matrix, fermionic_operators, generate_covariance

Z = np.array([[1, 0],
              [0, -1]], dtype='complex128')

def subsets(n,k):
    """returns all subsets of {1,...,n} of cardinality k"""
    
    return list(itertools.combinations(np.arange(1,n+1), k))
  
  
def ind(S):
    """given a list of integers >=1, S, returns a list with all entries shifted by -1"""
    
    return [(s-1) for s in S]


def HS(X,Y):
    """Hilbert-Schmidt inner product"""
    
    return np.trace(np.matmul(np.conj(X).T,Y))
 

def matching_sites(b,S):
    """given binary array b and index set S, outputs the number of entries of b with index in S that are 1"""
    
    if len(S) == 0:
        return 0
    
    else: 
        return np.sum(b[S])
    

def run_experiment(n, S, noise_channel, p, no_samples, no_trials):

	# define Gaussian state
	
    purestate = pure_initial(n)
    covar = generate_covariance(n)
    c, a = fermionic_operators(n)
    gaussianstate = gaussian_density_matrix(covar, c, a)
    
    
    # callibration procedure
    
    f_arr = []
    
    for k in range(0,n+1):
    	# fix callibration parameter f_2k
    	
    	# median of means
        estimates = np.zeros(no_trials)
        
        for trial in range(no_trials):
            print(trial)
            est = 0
            
            for sample in range(no_samples):
            
                Q, U = random_FGUclifford(n, ret_unitary=True)
                b = quantum_process(gaussianstate, U, n, p, False, noise_channel)
                est += get_f_2k(k,n,Q,b)

            est /= no_samples
            estimates[trial] = est
 
        f_arr.append(np.median(estimates))
	
	# estimation procedure
	
    estimates = np.zeros(no_trials)
    
    for trial in range(no_trials):
        print(trial)
        est = 0
        
        for sample in range(no_samples):
        
            Q, U = random_FGUclifford(n, ret_unitary=True)
            b = quantum_process(gaussianstate, U, n, p, False, noise_channel)
            est += estimate(n,f_arr,Q,b,S)
                
        est /= no_samples
        estimates[trial] = est
    
    return np.mean(estimates), f_arr




def get_f_2k(k,n,Q,b):
    """returns single-round estimator for f_2k, given measurement outcome b (an array) and matchgate Q"""
    
    sets = subsets(n,k)
    
    estimate = 0
    
    for S in sets:
        for Sp in sets:

            m = matching_sites(b, ind(Sp))
            
            estimate += 1/2**n*1/(scipy.special.binom(n,k))*(-1)**(k+m)*np.linalg.det(Q[np.ix_(ind(S), ind(Sp))])
            
    return estimate
    

def estimate(n,f_arr,Q,b,S):
    """returns single-round estimator for (-i)**(|S|/2)*gamma_S, given classical shadow (Q,b) and array of calibration parameters f_arr,
       It is not a problem to extend the function to compute expectation values of arbitrary observables, but this is a relevant special case (also computationally efficient)"""
   
    if len(S)%2==1:
    
        # parity preserved
        return 0
        
    k = len(S)//2
    sets = subsets(n,k)

    estimate = 0
    
    for Sp in sets:
        Sp_complete = np.array([(2*i-1,2*i) for i in Sp]).flatten()
         
        
        m = matching_sites(b, ind(Sp))

        estimate += 1/2**n*(-1)**(m)*np.linalg.det(Q[np.ix_(ind(S),ind(Sp_complete))])
        
    return (estimate)/f_arr[k] #(-1j)**k*(estimate)/f_arr[k] 
    
   
def true_val(n,S):
 
    purestate = pure_initial(n)
    covar = generate_covariance(n)
    c, a = fermionic_operators(n)
    gaussianstate = gaussian_density_matrix(covar, c, a)
    
    O = np.identity(2**n, dtype='complex128')
    
    for s in S:
        
        if (s-1)//2 == 0:
            
            if s ==1:
                op = X
                
            if s ==2:
                op = Y
            
            for qubit in range(1,n):
                op = np.kron(op, I)
        
        else:
            
            op = Z
            
            for qubit in range(1, s//2):
                op = np.kron(op, Z)
            
            if s%2==0:
                op = np.kron(op, Y)
            if s%2==1:
                op = np.kron(op, X)
            
            for qubit in range(s//2+1,n):
                op = np.kron(op, I)
            
            
        O = O @ op
        
    return HS(gaussianstate, O)


    
n = 3
S = [1,2]
noise_channel = "depolarizing"
p = 0
no_samples = 1000
no_trials = 2

#print(run_experiment(n, S, noise_channel, p, no_samples, no_trials))
true_val(n, S)



