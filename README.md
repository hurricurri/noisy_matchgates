# wrapper.py usage:

the main functions are "run_calibration" and 'run_estimation', which contains callibration/estimation procedures:

## run_calibration(n, noise_channel, p, no_samples_cali, no_trials_cali)
Input:
** n
(integer) the number of fermionic modes (or qubits, if one considers the Jordan-Wigner mapping)
*** noise_channel
(string) specifies which noise channel should be used
*** p
(float) the error probability
*** no_samples_cali and no_trials_cali
(integers) to obtain estimates for the callibration parameters $f_{2k}$, the median of means estimator is used.
this means that for every quantity that needs to be estimated, no_trials batches of no_samples single-shot estimates are obtained and the median of means is computed.

*** Return: f_arr
(array of floats) f_arr: callibration parameters (medians of means)

## run_estimation(n,p,noise_channel, f_arr, state, Q0,S0, no_samples_est, no_trials_est)

*** n
(integer) the number of fermionic modes (or qubits, if one considers the Jordan-Wigner mapping)
*** state
(numpy array) a 2^n times 2^n density matrix with respect to which one wants to estimate expectation values
*** S0
(ordered list of integers) specifies the observable $O$ for which one wants to estimate trace(gamma_S0 state)

example: for $S0=[1,2]$, $O$ would be $\gamma_1 \gamma_2$

All entries of $S0$ must be smaller than $2n$ (total number of majorana operators)
                  
*** noise_channel
(string) specifies which noise channel should be used

available are: 

"depolarizing"

"amplitude damping"

"X rotation"

*** p
(float) the error probability

example: for depolarizing noise, $p=1$ corresponds to fully depolarizing noise and $p=0$ corresponds to no noise.

*** no_samples_est and no_trials_est
(integers) to obtain estimates for $\operatorname{tr}(gamma_S0 state)$, the median of means estimator is used.
this means that for every quantity that needs to be estimated, no_trials batches of no_samples single-shot estimates are obtained and the median of means is computed.


*** outputs (complex number) median of means estimate for $\operatorname{tr}(O \rho)$

## in addition

the function "true_val" can be used to compute the true value of $\operatorname{tr}(O \rho)$. It takes as input the number of modes $n$, the state $\rho$,
and which product $\gamma_{S}$ one wants to consider. the list S is of the same form as in wrapper.py.


