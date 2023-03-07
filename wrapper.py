import numpy as np
import numpy.linalg as linalg
import itertools
import scipy.special
import matplotlib.pyplot as plt

from sample_matchgates import random_FGUclifford, I, X, Y
from quantum_process import quantum_process,noise_channel_ope
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

def run_calibration(n, noise_channel, p, no_samples, no_trials):
# callibration procedure
    print("callibration procedure")
    
    f_arr = []
    
    # initial state for callibration procedure
    all_zeros = np.zeros((2**n,2**n),dtype= "complex128")
    all_zeros[0,0] = 1
    
    for k in range(0,n+1):
        # fix callibration parameter f_2k
        
        # median of means
        estimates = np.zeros(no_trials)
        
        for trial in range(no_trials):
            #print('t',trial)
            est = 0
            
            for sample in range(no_samples):
            
                Q, U = random_FGUclifford(n, ret_unitary=True)
                
                b = quantum_process(all_zeros, U, n, p, False, noise_channel)
                
               # print('b',b)
                
                est += get_f_2k(k,n,Q,b)

            est /= no_samples
            estimates[trial] = est
         
        f_arr.append(np.median(estimates))
     #   print("f_" + str(2*k),'=',estimates)
    
    return f_arr

def run_estimation(n,p,noise_channel, f_arr, state, Q0,S0, no_samples, no_trials):
# estimation procedure
    """
    n: the number of qubits
    p: depolarizing noise strength
    """
    
    print("estimation procedure")
    
    estimates = np.zeros(no_trials, dtype='complex128')
    
    for trial in range(no_trials):
        #print('t',trial)
        est = 0
        
        for sample in range(no_samples):
        
            Q, U = random_FGUclifford(n, ret_unitary=True)
                
            b = quantum_process(state, U, n, p, False, noise_channel)
            
            est += estimate_gammaS(n,f_arr,Q0,Q,b,S0)
            #est += estimate_gammaS_U(n,f_arr,U0,U,b,S0)   
        est /= no_samples
        estimates[trial] = est
    
    return np.mean(estimates)
    
 """   
def run_experiment(n, state, S, noise_channel, p, no_samples, no_trials):    
    
    # callibration procedure
    print("callibration procedure")
    
    f_arr = []
    
    # initial state for callibration procedure
    all_zeros = np.zeros((2**n,2**n),dtype= "complex128")
    all_zeros[0,0] = 1
    
    for k in range(0,n+1):
        # fix callibration parameter f_2k
        
        # median of means
        estimates = np.zeros(no_trials)
        
        for trial in range(no_trials):
            #print('t',trial)
            est = 0
            
            for sample in range(no_samples):
            
                Q, U = random_FGUclifford(n, ret_unitary=True)
                
                b = quantum_process(all_zeros, U, n, p, False, noise_channel)
                
               # print('b',b)
                
                est += get_f_2k(k,n,Q,b)

            est /= no_samples
            estimates[trial] = est
         
        f_arr.append(np.median(estimates))
        print("f_" + str(2*k),'=',estimates)
    
    # estimation procedure
    print("estimation procedure")
    
    estimates = np.zeros(no_trials)
    
    for trial in range(no_trials):
        #print('t',trial)
        est = 0+0.j
        
        for sample in range(no_samples):
        
            Q, U = random_FGUclifford(n, ret_unitary=True)
           ### Q = np.identity(2 * n)
           # U = np.identity(2^n)
                
            b = quantum_process(state, U, n, p, False, noise_channel)
        
            est += estimate_gammaS(n,f_arr,Q0,Q,b,S0)
                
        est /= no_samples
        estimates[trial] = est
    
    return np.mean(estimates), f_arr
"""

def get_f_2k(k,n,Q,b):
    """returns single-round estimator for f_2k, given measurement outcome b (an array) and matchgate Q"""
    sets = subsets(n,k)
    estimate = 0
    for S in sets:
        temp=estimate
        S_complete = np.array([(2*i-1,2*i) for i in S]).flatten()
        for Sp in sets:
            Sp_complete = np.array([(2*i-1,2*i) for i in Sp]).flatten()
            m = matching_sites(b, ind(Sp))
            estimate += (-1)**(m)*np.linalg.det(Q[np.ix_(ind(Sp_complete), ind(S_complete))])
    estimate *= 1/(scipy.special.binom(n,k))
    #print(estimate)
    return estimate

"""
def get_f_2k_new(k,n,U,b):
    sets = subsets(n, k)
    est = 0
    all_zeros = np.zeros((2**n,2**n),dtype= "complex128")
    all_zeros[0,0] = 1
    num_b = 0
    temp = 1
    for j in range(n):
        num_b = num_b + temp * b[n-1-j]
        temp = temp*2
        
    mat_b = np.zeros((2**n,2**n),dtype= "complex128")
    mat_b[num_b,num_b] = 1 
    for S in sets:
        NS = np.array([(2*i-1,2*i) for i in S]).flatten()
        gamma_S = majorana_op(n,NS)
        gamma_S_new = np.matmul(U,np.conj(gamma_S).T)
        gamma_S_new = np.matmul(gamma_S_new, np.conj(U).T)
        temp = np.trace(np.matmul(all_zeros, gamma_S)) * np.trace(np.matmul(gamma_S_new, mat_b))
        est = est + temp
    est = est /(scipy.special.binom(n,k))
    
    return est"""

def estimate_gammaS_U(n,f_arr,U0,U,b,S0):
    """returns single-round estimator for tr(gamma_S rho) by U calculation"""
    k = len(S0)
    if k%2 == 1:
        return 0
    k=k//2
    gamma_S0 = majorana_op(n,S0)
    sets = subsets(2*n, 2*k)
    
    # array rep for b
    num_b = 0
    temp = 1
    for j in range(n):
        num_b = num_b + temp * b[n-1-j]
        temp = temp*2
    mat_b = np.zeros((2**n,2**n),dtype= "complex128")
    mat_b[num_b,num_b] = 1 
    
    res = 0
    for S in sets:
        gamma_S = majorana_op(n,S)
        temp= np.trace(np.conj(U0).T @ gamma_S0 @ U0 @ gamma_S)
        temp *= np.trace(U @ np.conj(gamma_S).T @ np.conj(U).T @ mat_b)
        res += temp
    res = res/f_arr[k]/2**n
    return res
    
def estimate_gammaS(n,f_arr,Q0,Q,b,S0):
    """returns single-round estimator for tr(gamma_S rho), given 
    the number of qubits n, 
    classical shadow (Q,b), 
    measurement gamma_S (Q0,S0)
    and array of calibration parameters f_arr"""
   
    if len(S0)%2==1:
    
        # parity preserved
        return 0
        
    k = len(S0)//2
    sets = subsets(n,k)
    sets0 = subsets(2*n,2*k)
    estimate = 0
    
    for S in sets0:
        
        temp_1 = np.linalg.det(Q0[np.ix_(ind(S0),ind(S))])
        if temp_1 == 0:
            continue
        temp_2 = 0
        for Sp in sets:
            Sp_complete = np.array([(2*i-1,2*i) for i in Sp]).flatten()
            m = matching_sites(b, ind(Sp))
            temp_2 += (-1)**(m)*np.linalg.det(Q[np.ix_(ind(Sp_complete),ind(S))])
           
        estimate += temp_1 * temp_2
   # print('b=',b,'Q=',Q,'temp_1=',temp_1,'temp_2=',temp_2)
   # print('before coef', estimate, 'after coef', 1.j**k*(estimate)/f_arr[k])
    estimate = 1.j**k*(estimate)/f_arr[k] 
    
    
    return estimate

def estimate_GaussianState(n,f_arr,Q,b,S):
    """returns single-round estimator for tr(rho_g rho), given classical shadow (Q,b) and array of calibration parameters f_arr,
    """
   
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
        
    return i**(k/2)*(estimate)/f_arr[k]

   
def true_val(n,state,S): #return true value
    
    O = majorana_op(n,S)
    
    return HS(state, O)


def majorana_op(n,S):#return gamma_S
    O = np.identity(2**n, dtype='complex128')
    
    if len(S) == 0:
        return O
     
    for s in S:
        
        if (s+1)//2 -1 == 0:
            
            if s ==1:
                op = X
                
            if s ==2:
                op = Y
            
            for qubit in range(1,n):
                op = np.kron(op, I)
        
        else:
            
            op = Z
            for qubit in range(1,(s+1)//2-1):
                op = np.kron(op, Z)
            
            if s%2==0:
                op = np.kron(op, Y)
            if s%2==1:
                op = np.kron(op, X)
            for qubit in range((s+1)//2,n):
                op = np.kron(op, I)
            
#        print('s',s,'op',op,'l_op',len(op),'l_O',len(O))
        O = O @ op
    
    return O

def to_binary(x,n): # convert x to n bit bitstring, e.g., (3,4)--> 0011
    temp = bin(x)[2:]
    x = ['0' for j in range(n)]
    start = n - len(temp)
    x[start:]=temp[:]
    for j in range(n):
        x[j] = int(x[j])
    return x

def double_set(S):
    DS = [0 for l in range(2*len(S))]
    k = 0
    for element in S:
        DS[2*k] = 2 * element-1
        DS[2*k+1] = 2 * element
        k = k + 1
    return DS

    
def ideal_f_2k_calculation(n,k, p, angle, noisechannel):
# return the ideal f_2k val
    S_set = subsets(n,k)
    l = len(S_set)
    val = 0
    for x in range(2**n):
        for S in S_set:
            x_arr = to_binary(x,n)
            temp = 0
 #           print('x_arr',x_arr)
            for cur_j in S:
                temp = temp + x_arr[int(cur_j-1)]
            temp = (-1)**temp
            rho_mat = [[0 for j in range(2**n)] for k in range(2**n)]
            rho_mat[x][x] = 1
            D_S = double_set(S)
 #           print('x',x,'S',S,'DS',D_S)
 #           print('rho',rho_mat,'O',majorana_op(n,D_S))
 #           print('sign',temp,'trace',np.trace(np.matmul(rho_mat,majorana_op(n,D_S))))
            noisy_ope = noise_channel_ope(majorana_op(n,D_S), n, p, angle, noisechannel)
            temp = temp * np.trace(np.matmul(rho_mat, noisy_ope))
            val = val + temp
    val = val * (-1.j)**k/(scipy.special.binom(2*n,2*k))/2**n
    return val

def Verify_gammaS_U(n,state,U0,S0):
    
    gamma_S0 = majorana_op(n,S0)
   # print(gamma_S0,'\nU0 dagger:',U0.conj().T,'transformed:',U0.conj().T @ gamma_S0 @ U0)
    gamma_after = np.matmul(np.conj(U0).T,gamma_S0)
    gamma_after = np.matmul(gamma_after,U0)
    gamma_after = np.matmul(gamma_after,state)
    res = np.trace(gamma_after)
    return res

def Verify_gammaS_Q(n,state,Q,S0):
    k = len(S0)
    S_set = subsets(2*n,k)
    res = 0
    for S in S_set:
        gamma_S = majorana_op(n,S)
        res += np.linalg.det(Q[np.ix_(ind(S0),ind(S))]) * np.trace(gamma_S @ state)
    return res
    
"""
n = 2
S0=[1,2]
#Q0, U0 = random_FGUclifford(n, ret_unitary=True)
U0 = np.identity(2**n, dtype='complex128')
Q0 = np.identity(2*n)
noise_channel = "depolarizing"
p = 0
###no_samples_arr = [100000]

no_samples_cali = 1000
no_trials_cali = 10

no_samples_est = 1000
no_trials_est = 10
est_arr = []

all_zeros = np.zeros((2**n,2**n),dtype= "complex128")
all_zeros[0,0] = 1"""


"""f_arr = run_calibration(n, noise_channel, p, no_samples_cali, no_trials_cali)
print('f_2k (0<=k<=n):',f_arr)
estimates = run_estimation(n,f_arr, all_zeros, Q0,S0, no_samples_est, no_trials_est)
est_ideal = Verify_gammaS(n,all_zeros,Q0,S0)
print('est=',estimates,'ideal=', est_ideal)"""


"""
for k in range(n+1):
    val = ideal_f_2k_calculation(n,k, p, False, noise_channel)
    print('n=',n,'f_', k,'=','val=',val)"""

#state = gaussian_state(n)
#print(run_experiment(n,state,S,noise_channel,p,no_samples,no_trials))


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

