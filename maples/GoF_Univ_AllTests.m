function [Sn,PV,Test] = GoF_Univ_AllTests(X,M)

% "Cramér-von Mises GoF tests for the Uniform (0,1), Normal, Exponential, 
%  Pareto and Rayleigh distributions"
% Input  -> X: n-dimensional data vetor
%           M: Number of parametric bootstrap samples 
% Output -> Sn: Vector of the five stats
%           PV: Estimated P-values of the five tests
%           Test: Vector of indicators ('1' for rejection, '0' otherwise)
% Necessary procedures: GoF_Univ_Uniform, GoF_Univ_Norm, GoF_Univ_Exp,
%                       GoF_Univ_Pareto, GoF_Univ_Rayleigh

Sn = zeros(1,5); PV = zeros(1,5); Test = zeros(1,5);

[Sn(1),PV(1),Test(1)] = GoF_Univ_Uniform(X,M);
[Sn(2),PV(2),Test(2)] = GoF_Univ_Norm(X,M);
[Sn(3),PV(3),Test(3)] = GoF_Univ_Exp(X,M);
[Sn(4),PV(4),Test(4)] = GoF_Univ_Pareto(X,M);
[Sn(5),PV(5),Test(5)] = GoF_Univ_Rayleigh(X,M);