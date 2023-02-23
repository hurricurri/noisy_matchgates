# wrapper.py usage:

the main function is "run_experiment", which combines callibration and estimation procedures 
(individual subtasks are outsourced to the other python scripts):

## inputs 
"run_experiment" takes the following inputs

### n
(integer) the number of fermionic modes (or qubits, if one considers the Jordan-Wigner mapping)
### state
(numpy array) a 2^n times 2^n density matrix with respect to which one wants to estimate expectation values
### S
(ordered list of integers) specifies the observable $O$ for which one wants to estimate $\operatorname{tr}(O \rho)$

example: for $S=[1,2]$, $O$ would be $\gamma_1 \gamma_2$

All entries of $S$ must be smaller than $2n$ (total number of majorana operators)
                  
### noise_channel
(string) specifies which noise channel should be used

available are: 

"depolarizing"

"amplitude damping"

"X rotation"

### p
(float) the error probability

example: for depolarizing noise, $p=1$ corresponds to fully depolarizing noise and $p=0$ corresponds to no noise.

### no_samples and no_trials
(integers) to obtain estimates for the callibration parameters $f_{2k}$ and $\operatorname{tr}(O \rho)$, the median of means estimator is used.
this means that for every quantity that needs to be estimated, no_trials batches of no_samples single-shot estimates are obtained and the median of means is computed.


## outputs

### (complex number) median of means estimate for $\operatorname{tr}(O \rho)$
    
### (array of floats) f_arr: callibration parameters (medians of means)

for $n$ modes, f_arr will contain $n+1$ callibration parameters (all $f_{2k}$ for $k$ between $0$ and $n$) 


## in addition

the function "true_val" can be used to compute the true value of $\operatorname{tr}(O \rho)$. It takes as input the number of modes $n$, the state $\rho$,
and which product $\gamma_{S}$ one wants to consider. the list S is of the same form as in wrapper.py.


