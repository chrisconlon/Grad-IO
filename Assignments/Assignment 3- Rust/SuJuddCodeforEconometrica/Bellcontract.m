function [EV, CbEV] = Bellcontract(thetaCost, TransProb, RC)

% This m-file solves the integrated Bellman equation using constraction
% mapping iteration in the NFXP algorithm. 
%
% source: Su and Judd (2011), Constrained Optimization Approaches to
% Estimation of Structural Models.
% Code Revised: Che-Lin Su, May 2010.

global x N M beta indices;

global EVold tol_inner BellEval;


EV0 = EVold;
            
% The cost function in Table X is: c(x) = 0.001*theta_1*x
Cost  = 0.001*thetaCost*x ; 

ii = 0;
norm = 1;

while norm > tol_inner 
    %  Let CbEV[i] represent - Cost[i] + beta*EV[i]; 
    %  this is the expected payoff at x[i] if regular maintenance is chosen
    
    CbEV = - Cost + beta*EV0;  
    CbEV(N+1:N+M)=CbEV(N);

    s1 = exp(CbEV(indices));
    s2 = exp(-RC+CbEV(1));
    s =  s1 + s2;
    logs = log(s);    
    EV = logs * TransProb;
    BellmanViola = abs(EV - EV0);
    norm = max(BellmanViola);
    EV0 = EV;
    ii = ii+ 1;
end

% BellEval is the number of contraction mapping iterations needed to solve
% the Bellman equation.
EVold = EV0;
BellEval = ii + BellEval;
