import numpy as np
import scipy.special
import itertools
from sample_matchgates import I, X, Y

Z = np.array([[1, 0],
              [0, -1]], dtype='complex128')

def subsets(n,k):
    """returns all subsets of {1,...,n} of cardinality k"""
    
    return list(itertools.combinations(np.arange(1,n+1), k))



def majorana_op(n,S):#return gamma_S
    O = np.identity(2**n, dtype='complex128')
    
    if S == [0]:
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

def to_binary(x,n): # convert x to a n bit bitstring, e.g., (3,4)--> 0011
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

    
def ideal_f_2k_calculation(n,k):
# return the ideal f_2k val
    S_set = subsets(n,k)
    l = len(S_set)
    val = 0
    for x in range(2**n):
        for S in S_set:
            x_arr = to_binary(x,n)
            temp = 0
            for cur_j in S:
                temp = temp + x_arr[int(cur_j-1)]
            temp = (-1)**temp
            rho_mat = [[0 for j in range(2**n)] for k in range(2**n)]
            rho_mat[x][x] = 1
            D_S = double_set(S)
 #           print('x',x,'S',S,'DS',D_S)
 #           print('rho',rho_mat,'O',majorana_op(n,D_S))
 #           print('sign',temp,'trace',np.trace(np.matmul(rho_mat,majorana_op(n,D_S))))
            temp = temp * np.trace(np.matmul(rho_mat,majorana_op(n,D_S)))
           
            val = val + temp
    val = val * (-1.j)**k/(scipy.special.binom(2*n,2*k))/2**n
    return val


n = 4
for k in range(n+1):
    val = ideal_f_2k_calculation(n,k)
    print('n=',n,'k=',k,'val=',val)