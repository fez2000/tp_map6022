function [Sn,PV,Test] = GoF_Univ_Rayleigh(X,M)

% "Cramér-von Mises GoF test for the Rayleigh distribution"
% Input  -> X: n-dimensional data vetor
%           M: Number of parametric bootstrap samples 
% Output -> Sn : Value of the test statistic
%           PV : Estimated P-value of the test
%           Test: '1' if Ho is rejected, '0' otherwise

% Parameter estimation
theta = (var(X) + mean(X)^2) / 2; 

% Computation of the test statistic
n = length(X);
SUM1 = 0; SUM2 = 0;
for i=1:n
    SUM2 = SUM2 + (0.5 - raylcdf(X(i),sqrt(theta))^2/2);
    for j=1:n
        SUM1 = SUM1 + (1 - raylcdf(max(X(i),X(j)),sqrt(theta)));
    end
end
Sn = (SUM1/n) - 2*SUM2 + (n/3);

% Estimation of the P-value
S = zeros(1,M);
parfor k=1:M
    SUM1 = 0; SUM2 = 0;
    R = raylrnd(sqrt(theta),1,n);
    theta_boot = (var(R) + mean(R)^2) / 2; 
    for i=1:n
        SUM2 = SUM2 + (0.5 - raylcdf(R(i),sqrt(theta_boot))^2/2);
        for j=1:n
            SUM1 = SUM1 + (1 - raylcdf(max(R(i),R(j)),sqrt(theta_boot)));
        end
    end
    S(k) = SUM1/n - 2*SUM2 + (n/3);
end
PV = sum(S>Sn) / M;
Test = (PV < 0.05);