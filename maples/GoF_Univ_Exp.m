function [Sn,PV,Test] = GoF_Univ_Exp(X,M)

% "Cramér-von Mises GoF test for the Exponential distribution"
% Input  -> X: n-dimensional data vetor
%           M: Number of parametric bootstrap samples 
% Output -> Sn : Value of the test statistic
%           PV : Estimated P-value of the test
%           Test: '1' if Ho is rejected, '0' otherwise

% Parameter estimation
lambda = mean(X);

% Computation of the test statistic
n = length(X);
SUM1 = 0; SUM2 = 0;  
for i=1:n
    SUM2 = SUM2 + (expcdf(X(i),lambda))^2;
    for j=1:n
        SUM1 = SUM1 + expcdf(max(X(i),X(j)),lambda);
    end
end
Sn = SUM2 - (SUM1/n) + (n/3);

% Estimation of the P-value
S = zeros(1,M);
parfor k=1:M
    SUM1 = 0; SUM2 = 0;
    R = exprnd(lambda,1,n);
    lambda_boot = mean(R); 
    for i=1:n
        SUM2 = SUM2 + (expcdf(R(i),lambda_boot))^2;
        for j=1:n
            SUM1 = SUM1 + expcdf(max(R(i),R(j)),lambda_boot);
        end
    end
    S(k) = SUM2 - (SUM1/n) + (n/3);
end
PV = sum(S>Sn) / M;
Test = (PV < 0.05);
   