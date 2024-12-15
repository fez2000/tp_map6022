function [Sn,PV,Test] = GoF_Univ_Pareto(X,M)

% "Cramér-von Mises GoF test for the Pareto distribution"
% Input  -> X: n-dimensional data vetor
%           M: Number of parametric bootstrap samples 
% Output -> Sn : Value of the test statistic
%           PV : Estimated P-value of the test
%           Test: '1' if Ho is rejected, '0' otherwise

% Parameter estimation
alpha = min(X); beta = (alpha/mean(X)) + 1;

% Computation of the test statistic
n = length(X);
SUM1 = 0; SUM2 = 0;  
for i=1:n
    SUM2 = SUM2 + ( 1 - (alpha/X(i))^beta )^2;
    for j=1:n
        SUM1 = SUM1 + 1 - (alpha/max(X(i),X(j)))^beta;
    end
end
Sn = SUM2 - (SUM1/n) + (n/3);

% Estimation of the P-value
S = zeros(1,M);
parfor k=1:M
    SUM1 = 0; SUM2 = 0;
    R = zeros(n,1);
    for i=1:n
        R(i) = alpha*(1-rand())^(-1/beta); 
    end
    alpha_boot = min(R); beta_boot = (alpha_boot/mean(R))+1; 
    for i=1:n
        SUM2 = SUM2 + ( 1 - (alpha_boot/R(i))^beta )^2;
        for j=1:n
            SUM1 = SUM1 + 1 - (alpha_boot/max(R(i),R(j)))^beta_boot;
        end
    end
    S(k) = SUM2 - (SUM1/n) + (n/3);
end
PV = sum(S>Sn) / M;
Test = (PV < 0.05);