function L = BivCopula_Likelihood(R,C,TauC)

% "Pseudo log-likelihood of a given copula"
% Input  -> R: n x 2 matrix of standardized ranks
%           C: From List#1b
%           TauC: Kendall's tau
% Output -> L: Value of the pseudo log-likelihood
% Necessary procedure: BivCopula_Density

n = length(R);
L = 0;
parfor i=1:n
    L = L + log(BivCopula_Density(C,TauC,R(i,1),R(i,2)));
end
L = real(L);