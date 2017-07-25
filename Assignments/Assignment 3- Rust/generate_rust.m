function [MC_dt,MC_xt,MC_dx]=simdata(param,EV)
%The following is the adjusted section of the code from Su and Judd RustBusMLETableX_MC.m

MC=param.MC;
beta = param.beta; 
nT = param.nT;
nBus = param.nBus;
N = param.N;
M = param.M;
RC = param.RC;
thetaCost = param.thetaCost;
thetaProbs = param.thetaProbs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Simulate Data for 250 date sets -- (state, decision) = (xt,dt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = (1:N)';
P0 = 1./ (1 + exp( 0.001*thetaCost*x - beta.*EV - RC - 0.001*thetaCost*x(1)+ beta*EV(1)));

rand_seed = 100;
rand('seed',rand_seed);

MC_xt = zeros(nT, nBus, MC);
MC_dt = zeros(nT, nBus, MC);
MC_dx = zeros(nT, nBus, MC);

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
    MC_dx(:,:,kk) = dx;

end