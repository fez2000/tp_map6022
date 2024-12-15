function [Sn,PV,Test] = GoF_Univ_Norm(X,M)

% "Cramér-von Mises GoF test for the Normal distribution"
% Input  -> X: n-dimensional data vetor
%           M: Number of parametric bootstrap samples 
% Output -> Sn : Value of the test statistic
%           PV : Estimated P-value of the test
%           Test: '1' if Ho is rejected, '0' otherwise

% Parameter estimation
mu = mean(X); sigma = std(X);

% Computation of the test statistic
n = length(X);
SUM1 = 0; SUM2 = 0;  
for i=1:n
    SUM2 = SUM2 + (normcdf(X(i),mu,sigma))^2;
    for j=1:n
        SUM1 = SUM1 + normcdf(max(X(i),X(j)),mu,sigma);
    end
end
Sn = SUM2 - (SUM1/n) + (n/3);

% Estimation of the P-value
S = zeros(1,M);
parfor k=1:M
    SUM1 = 0; SUM2 = 0;
    R = normrnd(mu,sigma,1,n);
    mu_boot = mean(R); sigma_boot = std(R); 
    for i=1:n
        SUM2 = SUM2 + (normcdf(R(i),mu_boot,sigma_boot))^2;
        for j=1:n
            SUM1 = SUM1 + normcdf(max(R(i),R(j)),mu_boot,sigma_boot);
        end
    end
    S(k) = SUM2 - (SUM1/n) + (n/3);
end
PV = sum(S>Sn) / M;
Test = (PV < 0.05);