function [ CbEV, EV, TransProb, thetaProbs, RC]  = DefVar(X)
% This file defines the decision variables for structural parameters and 
% expected value functions in the constrained optimization formulation 
% of the Harold Zurcher bus-engine replacement model.
% 
% source: Su and Judd (2011), Constrained Optimization Approaches to
% Estimation of Structural Models.
% Code Revised: Che-Lin Su, May 2010.

global x M beta;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DECLARE STRUCTURAL PARAMETERS TO BE ESTIMATED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define parameters for cost function
thetaCost = X(1);

% Define parameters and definition of transition process
% thetaProbs defines Markov chain
thetaProbs = X(2:M+1);
TransProb = thetaProbs;

% Define replacement cost parameter 
RC = X(7);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DECLARE EQUILIBRIUM CONSTRAINT VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The constrained optimization approach requires us to solve equilibrium 
% constraint variables. 
% Define expected value function variables
EV = X(8:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DECLARE AUXILIARY VARIABLES 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  Define auxiliary variables to economize on expressions	
%  Create Cost variable to represent the cost function; 
%  Cost[i] is the cost of regular maintenance at x[i].

% Cost function in Table X is: c(x) = 0.001*theta_1*x
Cost  = 0.001*thetaCost*x ; 

%  Let CbEV[i] represent - Cost[i] + beta*EV[i].  
%  This is the expected payoff at x[i] if regular maintenance is chosen
CbEV = - Cost + beta*EV;                     

end