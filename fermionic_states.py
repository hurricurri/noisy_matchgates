import numpy as np

##Create initial states

##Pure state with first orbital occupied
def pure_initial(n):
    rho = np.zeros((2**n, 1))
    rho[:1] = 1 
    density = rho @ rho.conj().T 
    return density


##Fermionic gaussian state 

def generate_covariance(n):
    '''Generates a covariance matrix for a fermionic wavefunction with n modes.'''

    #normal distribution
    correlations = np.random.normal(size=(n, n))

    #symmetrise the correlations matrix
    correlations = (correlations + correlations.T) / 2

    np.fill_diagonal(correlations, 0)

    #normalise
    covariance = np.dot(correlations, correlations.T) / (n - 1)

    return covariance


def fermionic_operators(n):

    creation = np.zeros(n, dtype='complex128')

    annihilation = np.zeros(n, dtype='complex128')

    for i in range(n-1):

        creation[i] = np.sqrt(i+1)

        annihilation[i+1] = np.sqrt(i+1)

    return creation, annihilation


def gaussian_density_matrix(covariance, creation, annihilation):

    n = creation.shape[0]

    eps = 1e-8 #small positive constant to ensure invertibility of covariance matrix

    covariance = covariance + eps * np.eye(n)

    a = np.dot(annihilation, np.linalg.cholesky(covariance))

    psi = np.exp(-0.5 * np.dot(a, a))

    rho = np.identity(2**n, dtype='complex128') * psi
    

    ##Jordan-Wigner mapping
    for i in range(n):

        X = np.identity(2**n, dtype='complex128')

        Z = np.identity(2**n, dtype='complex128')

        X[::2**(i+1), ::2**(i+1)] = 0

        X[2**(i)::2**(i+1), 2**(i)::2**(i+1)] = -1

        Z[::2**i, ::2**i] = -1

        for j in range(i):

            X[2**j::2**(i+1), 2**j::2**(i+1)] *= -1

        rho = (X @ rho @ X.conj().T) + (Z @ rho @ Z.conj().T)
        
        rho = rho/np.trace(rho) #normalize
        
    return rho

##Majorana operators state 
