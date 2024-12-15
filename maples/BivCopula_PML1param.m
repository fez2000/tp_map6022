function theta_HAT = BivCopula_PML1param(X,C)

% "Pseudo maximum likelihood (PML) estimator of a one-parameter copula"
% Input  -> X: n x 2 data matrix
%           C: Either C=1 (Normal), C=2 (centred Chi-square), C=3 (Student),
%              C=5 (Laplace p=1,2,3,4,5), C=9 (Plackett), C=11 (Clayton),
%              C=13 (Frank), C=15 (Gumbel) or C=17 (FGM)
% Output -> theta_HAT: PML estimator of the unknown parameter 
% Necessary procedures: VectorOfRanks, BivCopula_Likelihood,
%                       BivCopula_InversionKendall, fminsearchbnd

n = length(X); 
R = VectorOfRanks(X) / (n+1);

f = @(x)(-BivCopula_Likelihood(R,C,x));
LB = .01; UB = .99; x0 = (LB+UB)/2;
options = optimset('maxiter',25);
tau_HAT = fminsearchbnd(f,x0,LB,UB,options);

theta_HAT = BivCopula_InversionKendall(C,tau_HAT);
