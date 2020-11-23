# Load dependencies
import numpy as np
from numpy.random import default_rng

import scipy
import scipy.stats as stats
import scipy.sparse as sparse
from scipy.optimize import minimize
import pandas as pd

class RustModel(object):	
    # Constructor
    def __init__(self, x_min=0,x_max=200,n_grid=201,beta=0.975):
        self.n_grid = n_grid
        self.x_values = np.linspace(x_min,x_max,n_grid)/100.0

        # Parameters
        self.set_beta(beta)
        self.params = np.zeros(3)
        self.cost_param = np.zeros(2)
        self.RC = 0
        self.probs = np.array([1, 0, 0])

        # Parametrized Objects
        self.tpm = np.eye(n_grid)
        self.per_period_cost = np.zeros(n_grid)
        self.EV = np.zeros(n_grid)
        self.p_continue = np.zeros(n_grid)
        return 

    # Set cost parameters (quadratic plus replacement cost)
    def set_parameters(self,params=np.zeros(3)):
        if len(params) != 3:
            print("Parameter Vector Needs to be of length 3 (a x + b x ^2 , RC)")
        else:
            self.params = params
            self.cost_param = params[0:2]
            self.RC = params[2]
            self.update_cost_function()
        return

    # Set discount factor
    def set_beta(self,beta):
        if beta >= 1:
            print("Beta must be < 1!")
        self.beta = beta
        return 

    # Set the transition matrix
    def update_transition_probs(self,probs):
        # Start each row of transition probaiblities on the diagonal.
        A=sparse.diags(probs,range(0,len(probs)),shape=(self.n_grid,self.n_grid)).toarray()
        # Reminder: we go up to 4 states ahead of current state
        # which creates a problem towards boundary of state space and rows must sum to one.
        A[:,-1]+=1-A.sum(axis=1)
        self.tpm = A
        self.probs = np.array(probs)
        return

    def update_cost_function(self):
        # Quadratic function of the state variable 
        self.per_period_cost= self.x_values * self.cost_param[0] + (self.x_values**2) * self.cost_param[1]

    def solve_value_function(self):
        # Single Bellman Iteration
        # Note: we need a pure function of EV
        def fp_iter(EV,per_period_cost,RC):
            replace = self.beta * EV[0] -per_period_cost[0] - RC 
            wait    = self.beta * EV    -per_period_cost
            EV_new = self.tpm @ np.log(np.exp(replace)+np.exp(wait))
            return EV_new

        # Do value function iteration until we reach success
        # Warning: the accelerated version does del2 does something weird here
        EV_out=scipy.optimize.fixed_point(fp_iter,self.EV,maxiter=10000,method='iteration',args=(self.per_period_cost,self.RC))

        # Check that we succeeded
        if np.any(np.isnan(EV_out)):
            print("Value Function Iteration has failed.")
        else:
            self.EV=EV_out
        return

    def compute_choice_probabilities(self):
        # Given Ex-Ante Value function compute CCPs
        ex_ante_value = -self.per_period_cost + self.beta * self.EV
        value_diff = -ex_ante_value -self.RC + ex_ante_value[0]
        self.p_continue = 1./(1+np.exp(value_diff))
        
    # This code is needed only for generating fake data
    def simulate_bus(self,n_periods=120):
        n_probs = len(self.probs)
        rng=default_rng()
        u1=rng.random(n_periods) 
        
        # generate incremental mileage (+0,+1,+2,etc.)
        incremental_mileage=rng.choice(range(0,n_probs),size=n_periods,p=self.probs)

        # (x_it, y_it)
        state_var = np.zeros(n_periods).astype(int)
        replacement = np.zeros(n_periods).astype(int)
        counter = 0

        # increment mileage and check if we cross replacment threshold
        for idx, x in np.ndenumerate(incremental_mileage):
            counter += int(x)
            if self.p_continue[counter] < u1[idx]:
                replacement[idx]=1
                counter=x
            state_var[idx]=counter

        return np.vstack([range(0,n_periods),replacement,state_var]).transpose()