function theta_HAT = BivCopula_IKT(X,C)

% "Inversion of Kendall's tau (IKT) estimator of a given copula"
% Input  -> X: n x 2 data matrix
%           C: From List#1a
% Output -> theta_hat: IKT estimator
% Necessary procedures: KendallTau, BivCopula_InversionKendall

Tau_n = KendallTau(X);
theta_HAT = BivCopula_InversionKendall(C,Tau_n);
