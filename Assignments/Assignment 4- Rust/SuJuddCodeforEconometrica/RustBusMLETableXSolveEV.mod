# Source: Su and Judd (2011), Constrained Optimization Approaches to Estimation of Structural Models.
# Code Revised: Che-Lin Su, May 2010.
 
# The following AMPL files defines and solve the integrated Bellman equation for the expected value functions
# for given structural parameter values. 

# AMPL Model File:   RustBusMLETableXSolveEV.mod
# AMPL Data File:    RustBusMLETableXParamData.dat
# AMPL Command File: RustBusMLETableXSolveEV.run

# HAROLD ZURCHER BUS REPAIR EXAMPLE 
# The specifications of this model listed below follows those in Table X in Rust (1097, p.1022). 
#	Replacement Cost: RC
#       Cost Function c(x, theta1) = 0.001 * thetaCost * x
#       The mileage state space is discretized into 175 points (Fixed Point Dimension): x = 1, ... , 175 
#       Mileage transitions: move up at most t states. The transition probabilities are
# The data is simulated using the parameter values reported in Table X.

# SET UP THE MODEL and DATA #
  
#  Define and process the data
param nBus;               # number of buses in the data
set B := 1..nBus;         # B is the index set of buses
param nT;		  # number of periods in the data
set T := 1..nT;	          # T is the vector of time indices

#  Define the state space used in the dynamic programming part
param    N; 	          # number of states used in dynamic programming approximation
set X := 1..N; 	          # X is the index set of states
param x {i in X} := i;    #  x[i] denotes state i; 

# Parameters and definition of transition process
# In this example, M = 5: the bus mileage reading in the next period can either stay in current state or can move up to 4 states
param M;   

# Define discount factor. We fix beta since it can't be identified.
param beta;      	 # discount factor

# Data: (xt, dt) 
param dt {t in T, b in B};	 # decision of bus b at time t
param xt {t in T, b in B};       # mileage (state) of bus b at time t

# END OF MODEL and DATA SETUP #

# DEFINING STRUCTURAL PARAMETERS and ENDOGENOUS VARIABLES TO BE SOLVED #
# Parameters for (linear) cost function
#    c(x, thetaCost) = 0.001*thetaCost*x ;
var thetaCost >= 0; 

# thetaProbs[i] defines transition probability that mileage in next period moves up (i-1) grid point. M=5 in this example. 
var thetaProbs {1..M} >= 0; 

# Replacement cost
var RC >= 0;

# Define true structural parameter values
var RCsol;
var thetaCostsol;
var thetaProbssol {1..M};

# DECLARE EQUILIBRIUM CONSTRAINT VARIABLES 
# The NLP approach requires us to solve equilibrium constraint variables
var EV {X};        	# Expected Value Function of each state

# END OF DEFINING STRUCTURAL PARAMETERS AND ENDOGENOUS VARIABLES 

#  DECLARE AUXILIARY VARIABLES  #
#  Define auxiliary variables to economize on expressions	

#  Create Cost variable to represent the cost function; 
#  Cost[i] is the cost of regular maintenance at x[i].
var Cost {i in X} = 0.001*thetaCost*x[i];	    

#  Let CbEV[i] represent - Cost[i] + beta*EV[i]; 
#  this is the expected payoff at x[i] if regular maintenance is chosen
var CbEV {i in X} = - Cost[i] + beta*EV[i];    
                 
#  Let PayoffDiff[i] represent -CbEV[i] - RC + CbEV[1]; 
#  this is the difference in expected payoff at x[i] between engine replacement and regular maintenance
var PayoffDiff {i in X} = -CbEV[i] - RC + CbEV[1];               

#  Let ProbRegMaint[i] represent 1/(1+exp(PayoffDiff[i])); 
#  this is the probability of performing regular maintenance at state x[i];
var ProbRegMaint {i in X} = 1/(1+exp(PayoffDiff[i]));     

# BellmanViola represents violation of the Bellman equations. 
var BellmanViola {i in 1..(N-M+1)} = sum {j in 0..(M-1)} log(exp(CbEV[i+j])+ exp(-RC + CbEV[1]))* thetaProbs[j+1] - EV[i];

#  END OF DECLARING AUXILIARY VARIABLES #

# DEFINE OBJECTIVE FUNCTION AND CONSTRAINTS #

# Define objective function
# Since we are solving only for EV, we use 0 as the objective function.
maximize Likelihood0: 0 ;

#  Define the constraints

subject to 
#  Bellman equation for states below N-M
	Bellman_1toNminusM {i in X: i <= N-(M-1)}: 
		EV[i] = sum {j in 0..(M-1)} 
			log(exp(CbEV[i+j])+ exp(-RC + CbEV[1]))* thetaProbs[j+1];

#  Bellman equation for states above N-M, (we adjust transition probabilities to keep state in [xmin, xmax])
	Bellman_LastM {i in X: i > N-(M-1) and i <= N-1}: 
		EV[i] = (sum {j in 0..(N-i-1)} 
			log(exp(CbEV[i+j])+ exp(-RC + CbEV[1]))* thetaProbs[j+1]) 
			+ (1- sum {k in 0..(N-i-1)} thetaProbs[k+1]) * log(exp(CbEV[N])+ exp(-RC + CbEV[1]));

#  Bellman equation for state N
	Bellman_N:  EV[N] = log(exp(CbEV[N])+ exp(-RC + CbEV[1]));

#  The probability parameters in transition process must add to one
   # Probability: sum {i in 1..M} thetaProbs[i] = 1;

#  Put bound on EV; this should not bind, but is a cautionary step to help keep algorithm within bounds
    EVBound {i in X}: EV[i] <= 500;

# END OF DEFINING OBJECTIVE FUNCTION AND CONSTRAINTS

# DEFINE THE OPTIMIZATION PROBLEM #

# Name the problem
problem MPECZurcher:

# Choose the objective function
Likelihood0,


# List the variables
EV, RC, thetaCost, thetaProbs, Cost, CbEV, PayoffDiff, ProbRegMaint, BellmanViola, 


# List the constraints
Bellman_1toNminusM,
Bellman_LastM, 
Bellman_N, 
EVBound;

# END OF DEFINING THE MLE OPTIMIZATION PROBLEM
