function [Tn,PV] = GoF_EmpCopula(X,Co,N,B)

% "Cramér-von Mises GoF stat between the empirical copula and
%  a simulated version under a specified copula under Ho"
% Input  -> X: n x d data matrix
%           Co: Copula family under Ho
%               If d=2, from List#1a; if d>2, from List#2a
%           N: Number of samples for the simulated version
%           B: Number of parametric bootstrap samples 
% Output -> Sn : Test statistic
%           PV : Estimated P-value of the test
% Necessary procedures: VectorOfRanks, BivCopula_Simulation, KendallTau

[n,d] = size(X);

% Computation of the test stat
if d==2
    tauHAT = KendallTau(X);
    Xstar = BivCopula_Simulation(N,Co,tauHAT);
elseif d>2
    SigmaHAT = PairwiseCopula_IKT(X,Co);
    Xstar = PairwiseCopula_Simulation(N,Co,SigmaHAT);
end
R = VectorOfRanks(X) / n; Rstar = VectorOfRanks(Xstar) / N;

A1 = ones(n,n); A2 = ones(n,N); A3 = ones(N,N);
for j=1:d
   A1 = A1.*(1-max(R(:,j),R(:,j).'));
   A2 = A2.*(1-max(R(:,j),Rstar(:,j).'));
   A3 = A3.*(1-max(Rstar(:,j),Rstar(:,j).'));
end
SUM1 = sum(sum(A1)); SUM2 = sum(sum(A2)); SUM3 = sum(sum(A3));
Tn = (SUM1/n) - (2*SUM2/N) + (n*SUM3/N^2);

% Estimation of the P-value with parametric bootstrap samples
TnBoot = zeros(1,B);
if d==2
    parfor b=1:B
        Xb = BivCopula_Simulation(n,Co,tauHAT);
        tauBoot = KendallTau(Xb);
        Xbstar = BivCopula_Simulation(N,Co,tauBoot);
        R = VectorOfRanks(Xb) / n; Rstar = VectorOfRanks(Xbstar) / N;
    
        A1 = ones(n,n); A2 = ones(n,N); A3 = ones(N,N);
        for j=1:d
            A1 = A1.*(1-max(R(:,j),R(:,j).'));
            A2 = A2.*(1-max(R(:,j),Rstar(:,j).'));
            A3 = A3.*(1-max(Rstar(:,j),Rstar(:,j).'));
        end
        SUM1 = sum(sum(A1)); SUM2 = sum(sum(A2)); SUM3 = sum(sum(A3));
        TnBoot(b) = (SUM1/n) - (2*SUM2/N) + (n*SUM3/N^2);
    end
elseif d>2
    parfor b=1:B
        Xb = PairwiseCopula_Simulation(n,Co,SigmaHAT);
        SigmaBoot = PairwiseCopula_IKT(Xb,Co);
        Xbstar = PairwiseCopula_Simulation(N,Co,SigmaBoot);
        
        R = VectorOfRanks(Xb) / n; Rstar = VectorOfRanks(Xbstar) / N;
    
        A1 = ones(n,n); A2 = ones(n,N); A3 = ones(N,N);
        for j=1:d
            A1 = A1.*(1-max(R(:,j),R(:,j).'));
            A2 = A2.*(1-max(R(:,j),Rstar(:,j).'));
            A3 = A3.*(1-max(Rstar(:,j),Rstar(:,j).'));
        end
        SUM1 = sum(sum(A1)); SUM2 = sum(sum(A2)); SUM3 = sum(sum(A3));
        TnBoot(b) = (SUM1/n) - (2*SUM2/N) + (n*SUM3/N^2);
    end
end
PV = sum(TnBoot > Tn) / B;