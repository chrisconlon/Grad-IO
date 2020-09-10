function [f,g,h] = likelihood(X)

% This m-file computes the value and the gradient of the likelihood
% function in the constrained optimization formulation of the Harold 
% Zurcher bus-engine replacement model.
% f is the value of the likelihood function evaluated at the input variable X 
% g is the gradient of the likelihood function evaluated at the input variable X
% In this implementation, we do not supply second-order analytic
% derivatives. Hence, the hessian of the likelihood function, h, is empty [].
% 
% source: Su and Judd (2011), Constrained Optimization Approaches to
% Estimation of Structural Models.
% Code Revised: Che-Lin Su, May 2010. 

global dt xt nT PayoffDiffPrime TransProbPrime

[CbEV, EV, TransProb, thetaProbs, RC] = DefVar(X);

%  Let PayoffDiff[i] represent -CbEV[i] - RC + CbEV[1]; 
%  this is the difference in expected payoff at x[i] between engine replacement and regular maintenance
PayoffDiff  = -CbEV - RC + CbEV(1);               

%  Let ProbRegMaint[i] represent 1/(1+exp(PayoffDiff[i])); 
%  this is the probability of performing regular maintenance at state x[i];
ProbRegMaint = 1./(1+exp(PayoffDiff)); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  OBJECTIVE AND CONSTRAINT DEFINITIONS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define objective function: Likelihood function 

% The likelihood function contains two pieces
% First is the likelihood that the engine is replaced given time t state in the data.
% Second is the likelihood that the observed transition between t-1 and t
% would have occurred.

f = 0;                  % value of the likelihood function
g = zeros(length(X),1); % gradient of the likelihood funciton

xt2 = xt(2:nT,:);
xt1 = xt(1:nT-1,:);

dtMinus = (dt(2:nT,:)==0);
dtPlus = (dt(2:nT,:)==1);
dtM1Minus = (dt(1:nT-1,:)==0);
dtM1Plus = (dt(1:nT-1,:)==1);

% Constracut the value of the likelihood function
f1 = 1-ProbRegMaint(xt2(dtPlus));
f2 = ProbRegMaint(xt2(dtMinus));
f3 = TransProb( xt2( dtM1Plus ));
f4 = TransProb(xt2(dtM1Minus) - xt1(dtM1Minus)+1);

f = -( sum(log(f1))+ sum(log(f2))+ sum(log(f3))+ sum(log(f4))); 
       
% Construct the gradient of the likelihood function
d1 = PayoffDiffPrime(:,xt2(dtPlus))*(1-f1) ;        
d2 = -PayoffDiffPrime(:,xt2(dtMinus))*(1-f2);                   
d3 = TransProbPrime(:, xt2(dtM1Plus))*(1./f3);        
d4 = TransProbPrime(:, xt2(dtM1Minus)-xt1(dtM1Minus)+1)*(1./f4);  

g = -(d1+d2+d3+d4);   

if nargout > 2 
   h=[];
end
    
