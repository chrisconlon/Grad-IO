function [f, gNFXP, h] = likelihoodNFXP(theta)

% This m-file computes the value and the gradient of the likelihood
% function in the NFXP algorithm of the Harold Zurcher bus-engine replacement model.
% f is the value of the likelihood function evaluated at the structural parameter vector theta.  
% g is the gradient of the likelihood function evaluated at the structural parameter vector theta.
% In this implementation, we do not supply second-order analytic
% derivatives. Hence, the hessian of the likelihood function, h, is empty [].
% 
% source: Su and Judd (2011), Constrained Optimization Approaches to
% Estimation of Structural Models.
% Code Revised: Che-Lin Su, May 2010. 

global dt xt nT nBus N M M1 M2 PayoffDiffPrime TransProbPrime beta CbEVPrime indices;

% Define parameters for cost function
thetaCost = theta(1);

% Define parameters and definition of transition process
% thetaProbs defines Markov chain
thetaProbs = theta(2:6);
TransProb = thetaProbs;

% Define replacement cost parameter 
RC = theta(7);

ntheta = length(theta);

% Use constration mapping iteration to solve the integrated Bellman equations
[EV, CbEV] = Bellcontract(thetaCost, TransProb, RC);

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

f = 0;

g = zeros(length(theta)+N,1);
gPayoffPrime = g;
gTransProbPrime = g;

dtM1Minus = [];
dtM1Plus  = [];
dtMinus   = [];
dtPlus    = [];

for i = 1:nBus
    
    dtM1Minus = find(dt(1:(nT-1),i)==0);
    dtM1Plus  = find(dt(1:(nT-1),i)==1);
    dtMinus   = find(dt((2:nT),i)==0)+1;
    dtPlus    = find(dt((2:nT),i)==1)+1;

    ProbRegMaint(xt(dtPlus,i));
    ProbRegMaint(xt(dtMinus,i));
    TransProb( xt( dtM1Plus+1,i ) );
    TransProb(xt(dtM1Minus +1,i)-xt(dtM1Minus,i)+1);

    % Compute the value of the likelihood function
    f = f -( sum( log( 1-ProbRegMaint(xt(dtPlus,i))))...
             + sum( log( ProbRegMaint(xt(dtMinus,i)))) ...
             +  sum( log( TransProb( xt( dtM1Plus +1,i ) ) ))...
             + sum( log( TransProb(xt(dtM1Minus +1,i)-xt(dtM1Minus,i)+1))) ); 
   
    % Compute the gradient of the likelihood function     
    if nargout > 1
        d1 = PayoffDiffPrime(:,xt(dtPlus,i))*ProbRegMaint(xt(dtPlus,i)) ;
        d2 = - PayoffDiffPrime(:,xt(dtMinus,i))*( 1-ProbRegMaint(xt(dtMinus,i)) );
        d3 = TransProbPrime(:, xt( dtM1Plus +1,i ))*(1./TransProb( xt( dtM1Plus +1,i ) ));
        d4 = TransProbPrime(:, xt(dtM1Minus+1,i)-xt(dtM1Minus,i)+1)*(1./TransProb( xt(dtM1Minus +1,i)-xt(dtM1Minus,i)+1 ));
        
        gPayoffPrime = gPayoffPrime -(d1+d2);  
        gTransProbPrime = gTransProbPrime -(d3+d4);  
    end
    
end 

% Continue to compute the gradient of the likelihood function  
if nargout > 1
        
    gPayoffPrimetheta = gPayoffPrime(1:ntheta);
    gPayoffPrimeEV = gPayoffPrime(ntheta+1:end);
    gTransProbPrimetheta = gTransProbPrime(1:ntheta);
    
    s1 = exp(CbEV(indices));
    s2 = exp(-RC+CbEV(1));
    s =  s1 + s2;
    logs = log(s);
    
    Rprime = zeros(1,ntheta+N);
    Rprime(ntheta)=-1;
    Rprime(ntheta+1)=beta;
            
    d1 = ((CbEVPrime.*repmat(exp(CbEV),1,ntheta + N) + exp(-RC+CbEV(1))*repmat(Rprime,N+M,1)))./(repmat(exp(CbEV)+exp(-RC+CbEV(1)),1,ntheta + N));
       
    sum1 =  reshape(sum(reshape(d1(indices',:) .* repmat(repmat(TransProb,N,1),1,ntheta + N),M,N, ntheta + N )),N,ntheta + N); 
    sum2 = logs*TransProbPrime';
    TPrime = sum1 + sum2;
    
    EVPrime = [zeros(N,ntheta) eye(N)];
    dTdtheta = TPrime(:,1:ntheta);
    dTdEV = TPrime(:,(ntheta+1):ntheta+N);
    gNFXP = dTdtheta'*(inv(eye(N)-dTdEV))'*gPayoffPrimeEV + gPayoffPrimetheta + gTransProbPrimetheta;
end

if nargout > 2  
   h=[];
end