import numpy as np
import numpy.linalg as linalg
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

def gaussian_state(n):
    """output Gaussian state"""

    purestate = pure_initial(n)
    covar = generate_covariance(n)
    c, a = fermionic_operators(n)
    gaussianstate = gaussian_density_matrix(covar, c, a)
    
    return gaussianstate

def run_experiment(n, state, S, noise_channel, p, no_samples, no_trials):    
    
    # callibration procedure
    print("callibration procedure")
    
    f_arr = []
    
    # initial state for callibration procedure
    all_zeros = np.zeros((2**n,2**n),dtype= "complex128")
    all_zeros[0,0] = 1
    
    for k in range(0,n+1):
        # fix callibration parameter f_2k
        
        print("parameter " + str(2*k))
        
        # median of means
        estimates = np.zeros(no_trials)
        
        for trial in range(no_trials):
            print(trial)
            est = 0
            
            for sample in range(no_samples):
            
                Q, U = random_FGUclifford(n, ret_unitary=True)

                b = quantum_process(all_zeros, U, n, p, False, noise_channel)
                est += get_f_2k(k,n,Q,b)

            est /= no_samples
            estimates[trial] = est
 
        f_arr.append(np.median(estimates))
    
    # estimation procedure
    print("estimation procedure")
    
    estimates = np.zeros(no_trials, dtype="complex128")
    
    for trial in range(no_trials):
        print(trial)
        est = 0
        
        for sample in range(no_samples):
        
            Q, U = random_FGUclifford(n, ret_unitary=True)
            b = quantum_process(state, U, n, p, False, noise_channel)
            
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
            S_complete = np.array([(2*i-1,2*i) for i in S]).flatten()
            Sp_complete = np.array([(2*i-1,2*i) for i in Sp]).flatten()
        
            m = matching_sites(b, ind(Sp))
            
            estimate += 1/(scipy.special.binom(n,k))*(-1)**(m)*np.linalg.det(Q[np.ix_(ind(Sp_complete), ind(S_complete))])
    
    #print(estimate)
    
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

        estimate += (-1)**(m)*np.linalg.det(Q[np.ix_(ind(Sp_complete),ind(S))])
        
    return (+1j)**k*(estimate)/f_arr[k]
    
   
def true_val(n,state,S):
    
    O = majorana_op(n,S)
    
    return HS(state, O)


def majorana_op(n,S):
    O = np.identity(2**n, dtype='complex128')
    
    if S == [0]:
        return O
     
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
            
            for qubit in range(1, (s-1)//2):
                op = np.kron(op, Z)
            
            if s%2==0:
                op = np.kron(op, Y)
            if s%2==1:
                op = np.kron(op, X)
            
            for qubit in range((s-1)//2+1,n):
                op = np.kron(op, I)
                
        O = O @ op
    
    return O


n = 2
S = [1,2]
noise_channel = "depolarizing"
p = 0
#no_samples_arr = [100000]
no_samples = 10
#est_arr = []
no_trials = 15

state = gaussian_state(n)
print(run_experiment(n,state,S,noise_channel,p,no_samples,no_trials))
print(true_val(n, state, [1,2]))


#print(np.trace(state))
#print(np.allclose(state, np.conj(state).T, rtol=1e-05, atol=1e-08))
#print(linalg.eig(state)[0])

#print(majorana_op(2,[1]))
#print(majorana_op(2,[2]))
#print(majorana_op(2,[3]))
#print(majorana_op(2,[4]))

#Q, U = random_FGUclifford(n,ret_unitary = True)
#print(Q)
#print("#######################################")
#print(np.conj(U).T @ majorana_op(n,[0]) @ U)
#print("#######################################")
#print(np.conj(U).T @ majorana_op(n,[1]) @ U)
#print("#######################################")
#print(np.conj(U).T @ majorana_op(n,[2]) @ U)
#print("#######################################")
#print(np.conj(U).T @ majorana_op(n,[3]) @ U)
#print("#######################################")
#print(np.conj(U).T @ majorana_op(n,[4]) @ U)

#print(np.allclose(np.eye(U.shape[0]), np.matmul(np.conj(U).T, U)))

#for i, no_samples in enumerate(no_samples_arr):
#    print(run_experiment(n,state,S,noise_channel,p,no_samples,no_trials))
    #print("no_samples = " + str(no_samples))
    #est_arr.append(run_experiment(n,state,S,noise_channel,p,no_samples,no_trials)[0])

#print(true_val(n, state,S))
#plt.plot(no_samples_arr, est_arr, "o")
#plt.show()

#pot = 0
#k = 100000
#t = 2
#for i in range(k):
#    Q, U = random_FGUclifford(n,ret_unitary = True)
#    pot += 1/k*np.abs(np.trace(U))**(2*t)
#print(pot)

#print(majorana_op(3,[1,2,3,4]))

