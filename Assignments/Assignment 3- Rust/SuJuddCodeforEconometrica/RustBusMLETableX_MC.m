%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script file runs Monte Carlo experiments reported in Section 5 of
% the paper. 
%
% Users first need to specify the true values of structural parameters 
% as well as the discount factor in the data generating process. 
% Then the program will execute the following steps:
%
% 1) Call AMPL to solve the integrated Bellman equations for the expected 
%    value functions and compute the conditional choice probabilities to 
%    simulate data of mileage transitions and decisions for 250 data sets; 
%
% 2) Estimate the model using the constrained optimization approach with
%    AMPL/KNITRO implementation;
% 
% 3) Estimate the model using the constrained optimization approach with
%    MATLAB/ktrlink implementation with first-order analytic derivatives;
%
% 4) Estimate the model using the NFXP algorithm with MATLAB/ktrlink 
%    implementation with first-order analytic derivatives
%
% 5) Calculate summary statistics of the Monte Carlo experiments
%
% Source: Su and Judd (2011), Constrained Optimization Approaches to
% Estimation of Structural Models.
% Code Revised: Che-Lin Su, May 2010.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

diary RustBusMLETableX_MC_Multistart_beta975.out

global beta N M x xt dt nT nBus indices 
global PayoffDiffPrime TransProbPrime CbEVPrime RPrime EVPrime
global EVold tol_inner BellEval

% Specify the true values of structural parameters in the data generating
% process. For Monte Carlo experiments with a different the discount factor,
% change the value for beta below. 
beta = 0.975; 
nT = 120;
nBus = 50;
N = 175;
M = 5;
RC = 11.7257;
thetaCost = 2.4569;
thetaProbs = [ 0.0937
       0.4475
       0.4459
       0.0127
       0.0002 ];
   
thetatrue = [thetaCost; thetaProbs; RC];

MC = 250;   % number of monte carlo replications   
multistarts = 5; % number of starting points

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Step 1) 
%%% Call AMPL to solve the integrated Bellman equations for the expected 
%%% value functions and compute the conditional choice probabilities to 
%%% simulate data of mileage transitions and decisions for 250 data sets; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Call AMPL to solve for EV 

% Create the AMPL data file and write MATLAB data to AMPL format.
fid = fopen('RustBusMLETableXSolveEV.dat', 'w');    
    fprintf (fid, '# Data generated: %s\n\n', datestr(now));  
    fprintf (fid, 'data; \n\n'); 
    fprintAmplParamCLSU(fid, 'N', N, 1);
    fprintAmplParamCLSU(fid, 'M', M, 1);
    fprintAmplParamCLSU(fid, 'nBus', nBus, 1);
    fprintAmplParamCLSU(fid, 'nT', nT, 1);
    fprintAmplParamCLSU(fid, 'beta', beta, 1);
    fprintAmplParamCLSU(fid, 'RCsol', RC, 1);
    fprintAmplParamCLSU(fid, 'thetaCostsol',thetaCost, 1);
    fprintAmplParamCLSU(fid, 'thetaProbssol', thetaProbs, 1);
fclose(fid);
% The ampl executable in my computer is in the directory:
% "/Applications/AMPL/ampl64"
% Change the path below to the location of your AMPL executable.
strAmplCommand = '/Applications/AMPL/ampl64/ampl';    
outname = ['EV/SolveEV.out'];  
strAmplSystemCall = sprintf('%s RustBusMLETableXSolveEV.run > %s', strAmplCommand, outname);
[status,result] = system(strAmplSystemCall);
EV = csvread('EV/EV.sol');

figure;
plot(EV);
title(['Ev(x,0) for \beta =' num2str(beta) '; Fixed Point Dim = 175']);
xlabel('x');
ylabel('EV(x,0)');

save (['truethetaEV_beta' num2str(1000*beta)], 'beta', 'nT', 'nBus', 'N', 'M', 'RC', 'thetaCost', 'thetaProbs', 'EV', 'x');   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Simulate Data for 250 date sets -- (state, decision) = (xt,dt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = (1:N)';
P0 = 1./ (1 + exp( 0.001*thetaCost*x - beta.*EV - RC - 0.001*thetaCost*x(1)+ beta*EV(1)));

rand_seed = 100;
rand('seed',rand_seed);

MC_xt = zeros(nT, nBus, MC);
MC_dt = zeros(nT, nBus, MC);

for kk = 1:MC
   
    Rx  = unifrnd(0, 1, nT, nBus);
    Rd  = unifrnd(0, 1, nT, nBus);

    xt = ones(nT, nBus);
    dx = zeros(nT, nBus);
    dt = zeros(nT, nBus);
    cumTransProb = cumsum(thetaProbs);

    for t = 1:nT
        dt(t,:) = (Rd(t,:) >= P0(xt(t,:))');
        for i = 1:nBus
            dx(t,i) = find(Rx(t,i) < cumTransProb,1);
            if t < nT
                if dt(t,i) == 1
                   xt(t+1,i) = 1 + dx(t,i)-1;
                else 
                   xt(t+1,i) = min(xt(t,i) + dx(t,i)-1,N);
                end
            end
        end
    end
    MC_dt(:,:,kk) = dt;
    MC_xt(:,:,kk) = xt;
end

save (['RustBusTableXSimDataMC' num2str(MC) '_beta' num2str(1000*beta)], 'beta', 'nT', 'nBus', 'N', 'M', 'RC', 'thetaCost', 'thetaProbs', 'EV', 'x', 'MC', 'MC_dt', 'MC_xt');   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Setup Optimization Problem for MPEC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Construct the sparsity pattern of the constraint Jacobian and Hessian
        
PayoffDiffPrime = zeros(7+N,length(x));        
PayoffDiffPrime(1,:)=0.001*(x'-repmat(x(1),1,length(x)));       
PayoffDiffPrime(7,:)=-1;        
PayoffDiffPrime(8,:)= beta;
                
PayoffDiffPrime= -beta*[zeros(7,N); eye(N)] + PayoffDiffPrime;
                
CbEVPrime = zeros(length(x),7+N);        
CbEVPrime(:,1)=-0.001*x;        
                
CbEVPrime= beta*[zeros(N,7) eye(N)] + CbEVPrime;        
CbEVPrime =  [ CbEVPrime; repmat(CbEVPrime(length(x),:),M,1)];
                
TransProbPrime = zeros(7+N,M);        
TransProbPrime(2:6,:) = eye(M);           

RPrime = zeros(1,8-1+N);
RPrime(:,7)=-1;
RPrime(:,8)=beta;
    
indices = repmat((1:N)',1,M)+repmat((1:M),N,1)-1;
d1 = ((CbEVPrime.*repmat(ones(N+M,1),1,N+7) + repmat(RPrime,N+M,1)))./(repmat(ones(N+M,1),1,N+7));
       
sum1 =  reshape(sum(reshape(d1(indices',:) .* repmat(repmat(ones(M,1),N,1),1,N+7),M,N, N+7 )),N, N+7); 
sum2 = ones(N,M)*TransProbPrime';
    
EVPrime = [zeros(N,7) eye(N)];
    
JacobSpaPattern = (sum1 + sum2 - EVPrime);
JacobSpaPattern = ceil(abs(JacobSpaPattern))./max(ceil(abs(JacobSpaPattern)),1);
    
HessSpaPattern = ones(8,7+N);
HessSpaPattern = [HessSpaPattern; ones(N-1,8) eye(N-1)];
         
%Define upper and lower bounds for the following decision varialbes: 
%thetaCost,thetaProbs, and EV

%thetaCost
lb = zeros(N+7,1);
ub = zeros(N+7,1);

lb(1)=0;
ub(1)=inf;

%thetaProbs
lb(2:6)=0;
ub(2:6)=1;

%RC
lb(7)=0;
ub(7)=inf;

%EV
%Put bound on EV; this should not bind, but is a cautionary step to help keep algorithm within bounds
%lb(7:(7-1+N))=0;
lb(8:end)=-inf;
ub(8:end)=50;

%The probability parameters in transition process must add to one
%Linear constraint : sum thetaProbs = 1
Aeq = zeros(1,N+7);
Aeq(2:6)=1;
beq=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Setup Optimization Problem for NFXP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Define upper and lower bounds for structural parameters:
% thetaCost, thetaProbs

%thetaCost
NFXPlb = zeros(7,1);
NFXPub = zeros(7,1);

NFXPlb(1)=0;
NFXPub(1)=inf;

%thetaProbs
NFXPlb(2:6)=0;
NFXPub(2:6)=1;
NFXPlb(7)=0;
NFXPub(7)=inf;

%The probability parameters in transition process must add to one
%Linear constraint : sum thetaProbs = 1
NFXPAeq = zeros(1, length(thetatrue));
NFXPAeq(2:6)=1;
NFXPbeq=1;

%----------------------

KnitroExitAMPL = -10000*ones(MC,1);
thetaCostAMPL = zeros(1,MC);
RCAMPL = zeros(1, MC);
EVAMPL = zeros(N, MC);
thetaProbsAMPL = zeros(M,MC);
SolveTimeAMPL = zeros(MC,1);
ObjValAMPL = zeros(MC,1);
IterAMPL = zeros(MC,1);
FunEvalAMPL = zeros(MC,1);
SuccessAMPL = zeros(MC,1);

MC_status = [];
MC_result = [];

thetaMPECsol = zeros(length(thetatrue),MC);
tMPECsol = zeros(MC,1);
fvalMPECsol = zeros(MC,1);
flagMPECsol = -10000*ones(MC,1);
IterMPECsol = zeros(MC,1);
FunEvalMPECsol = zeros(MC,1);
SuccessMPEC = zeros(MC,1);

thetaNFXPsol = zeros(length(thetatrue),MC);
tNFXPsol = zeros(MC,1);
fvalNFXPsol = zeros(MC,1);
flagNFXPsol = -100000*ones(MC,1);
IterNFXPsol = zeros(MC,1);
FunEvalNFXPsol = zeros(MC,1);
SuccessNFXP = zeros(MC,1);
numBellEvalsol = zeros(MC,1);


%Generating starting Points
X0 = zeros(7+N,multistarts);
X0(1,:)= (1:1:multistarts);
X0(2:6,:)=1/M;
X0(7,:)= (4:1:4+(multistarts-1));

%Initial estimates for NFXP
% theta0 = X0(1:7);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Step 2) 
%%% Estimate the model using the constrained optimization approach with
%%% AMPL/KNITRO implementation;
%%%
%%% IMPLEMENTATION 1: MPEC/AMPL
%%% We implement the MPEC approach using AMPL modeling language
%%% with KNITRO as the solver.
%%% AMPL supplies first-order and second-order analytical derivatives
%%% and the sparsity pattern of the constraint Jacobian and the Hessian
%%% to the solver. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for kk = 1:MC
   
    fprintf('This is Monte Carlo run #%d out of %d\n \n', kk, MC);
    
    xt = MC_xt(:,:,kk);
    dt = MC_dt(:,:,kk);
    
    G = -1.0e10;
    
    for reps = 1:multistarts 
        
        fprintf('Starting MPEC/AMPL/KNITRO in Monte Carlo run #%d out of %d replications ... \n', kk, MC);
        disp(['Running Starting Point #' num2str(reps)])
        
        fid = fopen('RustBusMLETableX.dat', 'w');
        fprintf (fid, '# Data generated: %s\n\n', datestr(now));  
        fprintf (fid, 'data; \n\n'); 
        fprintAmplParamCLSU(fid, 'N', N, 1);
        fprintAmplParamCLSU(fid, 'M', M, 1);
        fprintAmplParamCLSU(fid, 'nBus', nBus, 1);
        fprintAmplParamCLSU(fid, 'nT', nT, 1);
        fprintAmplParamCLSU(fid, 'beta', beta, 1);
        fprintAmplParamCLSU(fid, 'inithetaCost', X0(1,reps), 1);
        fprintAmplParamCLSU(fid, 'iniRC', X0(7,reps), 1);
        fprintAmplParamCLSU(fid, 'iniEV', X0(8,reps), 1);
        fprintAmplParamCLSU(fid, 'xt', xt, 1); 
        fprintAmplParamCLSU(fid, 'dt', dt, 1);
        fclose(fid);
        % The ampl executable in my computer is in the directory:
        % "/Applications/AMPL/ampl64"
        % Change the path below to the location of your AMPL executable.
        strAmplCommand = '/Applications/AMPL/ampl64/ampl'; 
        outname = ['output/MC' num2str(kk) '_multistart' num2str(reps) '.out']; 
        strAmplSystemCall = sprintf('%s RustBusMLETableX.run > %s', strAmplCommand, outname);
        [status,result] = system(strAmplSystemCall);
    
        MC_status = [MC_status; status];
        MC_result = [MC_result; result];
        
        SolveTimeAMPL_reps = csvread('output/solvetime.sol');
        SolveTimeAMPL(kk) = SolveTimeAMPL(kk) + SolveTimeAMPL_reps;
        
        fid = fopen('KnitroMessage.sol');
        line1 = fgetl(fid);
        line2 = fgetl(fid);
        line3 = fgetl(fid);
    
        IterAMPL_reps = sscanf(line3, '%d',1);
        FunEvalAMPL_reps = sscanf(line3, '%*s %*s %d',1);
        fclose(fid); 
        
        fprintf('\n');
        disp([line1])
        disp([line2])
        disp([line3])
        fprintf('\n');
        
        IterAMPL(kk) = IterAMPL(kk) + IterAMPL_reps;
        FunEvalAMPL(kk) = FunEvalAMPL(kk) + FunEvalAMPL_reps;
        
        ObjVal_reps = csvread('output/objval.sol');
        KnitroExitFlag_reps = csvread('output/KnitroExit.sol'); 
        
        if KnitroExitFlag_reps == 0
            SuccessAMPL(kk) = SuccessAMPL(kk)+1;
            KnitroExitAMPL(kk) = KnitroExitFlag_reps;
            if ObjVal_reps > G   
                ObjValAMPL(kk) = ObjVal_reps;
                thetaCostAMPL(kk) = csvread('output/thetaCost.sol');
                thetaProbsAMPL(:,kk) = csvread('output/thetaProbs.sol');
                RCAMPL(kk) = csvread('output/RC.sol');
                EVAMPL(:,kk) = csvread('output/EV.sol');
                G = ObjVal_reps;
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Step 3) 
    %%% Estimate the model using the constrained optimization approach with
    %%%  MATLAB/ktrlink implementation with first-order analytic derivatives;
    %%%
    %%% IMPLEMENTATION 2: MPEC/MATLAB
    %%% We implement the MPEC approach using MATLAB programming language
    %%% with KNITRO (ktrlink) as the solver. 
    %%% We provide first-order analytical derivatives and sparsity pattern
    %%% of the constraint Jacobian to the solver.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ktroptsMPEC = optimset('DerivativeCheck','on','Display','iter',...
        'GradConstr','on','GradObj','on','TolCon',1E-6,'TolFun',1E-6,'TolX',1E-15,'JacobPattern',JacobSpaPattern, 'HessPattern', HessSpaPattern); 

    G = -1.0e10;
    
    for reps = 1:multistarts
        fprintf('Starting MPEC/ktrlink in Monte Carlo run #%d out of %d replications... \n', kk, MC);
        disp(['Running Starting Point #' num2str(reps)]) 
        x0 = X0(:,reps);
        t1 = cputime;    
        % profile on;
        [XMPEC_reps fvalMPEC_reps flagMPEC_reps outputMPEC] = ktrlink(@likelihood,x0,[],[],Aeq,beq,lb,ub,@confuneq,ktroptsMPEC,'knitroOptions.opt');  
        % profile viewer;  
        tMPEC_reps = cputime - t1;
        
        tMPECsol(kk) = tMPECsol(kk) + tMPEC_reps;
        IterMPECsol(kk) = IterMPECsol(kk) + outputMPEC.iterations;
        FunEvalMPECsol(kk) = FunEvalMPECsol(kk) + outputMPEC.funcCount;
        
        if flagMPEC_reps == 0
            SuccessMPEC(kk) = SuccessMPEC(kk)+1;
            flagMPECsol(kk) = flagMPEC_reps;        
            if -fvalMPEC_reps > G
                fvalMPECsol(kk) = -fvalMPEC_reps;
                thetaMPECsol(:,kk) = XMPEC_reps(1:7);
                G = -fvalMPEC_reps;
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Step 4) 
    %%% Estimate the model using the NFXP algorithm with MATLAB/ktrlink 
    %%% implementation with first-order analytic derivatives
    %%%
    %%% IMPLEMENTATION 3: NFXP/MATLAB
    %%% We implement the NFXP algorithm using MATLAB programming language
    %%% with KNITRO (ktrlink) as the solver. 
    %%% We provide first-order analytical derivatives to the solver.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ktroptsNFXP = optimset('Display','iter','GradObj','on','TolCon',1E-6,'TolFun',1E-6,'TolX',1E-15);
    
    G = -1.0e10;
    
    for reps = 1:multistarts
        fprintf('Running NFXP in Monte Carlo run #%d out of %d replications... \n', kk, MC); 
        disp(['Running Starting Point #' num2str(reps)]) 
        EVold = zeros(N,1);
        tol_inner = 1.e-10;
        BellEval = 0;
        theta0 = X0(1:7,reps);   
        t2 = cputime;  
        [thetaNFXP_reps fvalNFXP_reps flagNFXP_reps outputNFXP] = ktrlink(@likelihoodNFXP,theta0,[],[],NFXPAeq,NFXPbeq,NFXPlb,NFXPub,[],ktroptsNFXP,'knitroOptions.opt');     
        tNFXP_reps = cputime - t2;
        tNFXPsol(kk) = tNFXPsol(kk) + tNFXP_reps;
        IterNFXPsol(kk) = IterNFXPsol(kk) + outputNFXP.iterations;
        FunEvalNFXPsol(kk) = FunEvalNFXPsol(kk) + outputNFXP.funcCount;
        numBellEvalsol(kk) = numBellEvalsol(kk) + BellEval;
        
        if flagNFXP_reps == 0
            SuccessNFXP(kk) = SuccessNFXP(kk)+1;
            flagNFXPsol(kk) = flagNFXP_reps;
            if -fvalNFXP_reps > G
                fvalNFXPsol(kk) = -fvalNFXP_reps;
                thetaNFXPsol(:,kk) = thetaNFXP_reps;
                G = -fvalNFXP_reps;
            end
        end
    end
end


thetaAMPLsol = [thetaCostAMPL; thetaProbsAMPL; RCAMPL];

save(['MC' num2str(MC) '_beta' num2str(1000*beta) '_result'], 'thetatrue', 'thetaAMPLsol', 'thetaMPECsol', 'thetaNFXPsol', 'KnitroExitAMPL', ...
    'SolveTimeAMPL', 'ObjValAMPL', 'IterAMPL', 'FunEvalAMPL', 'SuccessAMPL','flagMPECsol','tMPECsol', 'fvalMPECsol', 'IterMPECsol', ...
    'FunEvalMPECsol', 'SuccessMPEC', 'flagNFXPsol','tNFXPsol', 'fvalNFXPsol', 'IterNFXPsol', 'FunEvalNFXPsol', 'numBellEvalsol', 'SuccessNFXP');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Step 5
%%% Calculate summary statistics of the Monte Carlo experiments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TotalSuccessAMPL = sum(SuccessAMPL);
fprintf('MPEC/AMPL finishes successfully in #%d out of %d runs \n \n', TotalSuccessAMPL, MC*multistarts);
meanthetaAMPL = mean(thetaAMPLsol(:, KnitroExitAMPL==0),2);
stdevthetaAMPL = std(thetaAMPLsol(:, KnitroExitAMPL==0),1,2);
biasthetaAMPL = meanthetaAMPL - thetatrue;
RMSEthetaAMPL = sqrt(mean((thetaAMPLsol(:, KnitroExitAMPL==0)- repmat(thetatrue,1,sum(KnitroExitAMPL==0))).^2,2));
meanObjValAMPL = mean(ObjValAMPL(KnitroExitAMPL==0));

meantAMPL = mean(SolveTimeAMPL)/multistarts;
meanIterAMPL = mean(IterAMPL)/multistarts;
meanFunEvalAMPL = mean(FunEvalAMPL)/multistarts;

TotalSuccessMPEC = sum(SuccessMPEC);
fprintf('MPEC/ktrlink finishes successfully in #%d out of %d runs \n \n', TotalSuccessMPEC, MC*multistarts);

meanthetaMPEC = mean(thetaMPECsol(:,flagMPECsol==0),2);
stdevthetaMPEC = std(thetaMPECsol(:,flagMPECsol==0),1,2);
biasthetaMPEC = meanthetaMPEC - thetatrue;
RMSEthetaMPEC = sqrt(mean((thetaMPECsol(:, flagMPECsol==0) - repmat(thetatrue,1,sum(flagMPECsol==0))).^2,2));
meanObjValMPEC = mean(fvalMPECsol(flagMPECsol==0));

meantMPEC = mean(tMPECsol)/multistarts;
meanIterMPEC = mean(IterMPECsol)/multistarts; 
meanFunEvalMPEC = mean(FunEvalMPECsol)/multistarts;

TotalSuccessNFXP = sum(SuccessNFXP);
fprintf('NFXP/ktrlink finishes successfully in #%d out of %d runs \n \n', TotalSuccessNFXP, MC*multistarts);

meanthetaNFXP = mean(thetaNFXPsol(:,flagNFXPsol==0),2);
stdevthetaNFXP = std(thetaNFXPsol(:,flagNFXPsol==0),1,2);
biasthetaNFXP = meanthetaNFXP - thetatrue;
RMSEthetaNFXP = sqrt(mean((thetaNFXPsol(:, flagNFXPsol==0)- repmat(thetatrue,1,sum(flagNFXPsol==0))).^2,2));
meanObjValNFXP = mean(fvalNFXPsol(flagNFXPsol==0));

meantNFXP = mean(tNFXPsol)/multistarts;
meanIterNFXP = mean(IterNFXPsol)/multistarts; 
meanFunEvalNFXP = mean(FunEvalNFXPsol)/multistarts;
meanBellEvalNFXP = mean(numBellEvalsol)/multistarts;

save(['MC' num2str(MC) '_beta' num2str(1000*beta) '_summary'], 'thetatrue', 'TotalSuccessAMPL', 'meanthetaAMPL', 'stdevthetaAMPL', 'biasthetaAMPL', 'RMSEthetaAMPL', 'meanObjValAMPL','meantAMPL', 'meanIterAMPL','meanFunEvalAMPL', ...
    'TotalSuccessMPEC', 'meanthetaMPEC', 'stdevthetaMPEC', 'biasthetaMPEC', 'RMSEthetaMPEC', 'meanObjValMPEC','meantMPEC', 'meanIterMPEC','meanFunEvalMPEC', ...
    'TotalSuccessNFXP', 'meanthetaNFXP', 'stdevthetaNFXP', 'biasthetaNFXP', 'RMSEthetaNFXP', 'meanObjValNFXP','meantNFXP', 'meanIterNFXP','meanFunEvalNFXP', 'meanBellEvalNFXP');

fprintf('The truth, mean estimates of MPEC/AMPL, mean estimates of MPEC/ktrlink and mean estimates of NFXP/ktrlink are: \n');
[thetatrue meanthetaAMPL meanthetaMPEC meanthetaNFXP] 

fprintf('The std deviation of estimates of MPEC/AMPL, MPEC/ktrlink, and NFXP/ktrlink are: \n');
[stdevthetaAMPL stdevthetaMPEC stdevthetaNFXP] 

fprintf('The bias in estimates of MPEC/AMPL, MPEC/ktrlink, and NFXP/ktrlink are: \n');
[biasthetaAMPL biasthetaMPEC biasthetaNFXP] 

fprintf('The RMSE of estimates of MPEC/AMPL, MPEC/ktrlink, and NFXP/ktrlink are: \n');
[RMSEthetaAMPL RMSEthetaMPEC RMSEthetaNFXP] 

fprintf('The average computational time (in seconds) in each run for MPEC/AMPL, MPEC/ktrlink, and NFXP/ktrlink are: \n');
[meantAMPL meantMPEC meantNFXP] 

fprintf('The average # of iterations in each run for MPEC/AMPL, MPEC/ktrlink, and NFXP/ktrlink are: \n');
[meanIterAMPL meanIterMPEC meanIterNFXP] 

fprintf('The average # of function evaluations in each run for MPEC/AMPL, MPEC/ktrlink, and NFXP/ktrlink are: \n');
[meanFunEvalAMPL meanFunEvalMPEC meanFunEvalNFXP] 

fprintf('The average # of Bellman iterations in the inner loop of NFXP in each run of NFXP/ktrlink is: \n');
[meanBellEvalNFXP] 

diary off;
