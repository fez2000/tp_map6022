function rho_n = SpearmanRho(X)

% "Sample value of Spearman's rank correlation"
% Input  -> X: n x 2 data matrix
% Output -> rho_n: Empirical Spearman's rho
% Necessary procedure: VectorOfRanks

n = length(X); 
R = VectorOfRanks(X);
SUM = (R(:,1)-R(:,2)).'*(R(:,1)-R(:,2));
rho_n = 1 - (6/(n*(n^2-1)))*SUM;
