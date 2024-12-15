function [Sn,PV,Test] = GoF_Univ_Uniform(U,M)

% "Cramér-von Mises GoF test for the Uniform distribution"
% Input  -> X: n-dimensional data vetor
%           M: Number of parametric bootstrap samples 
% Output -> Sn : Value of the test statistic
%           PV : Estimated P-value of the test
%           Test: '1' if Ho is rejected, '0' otherwise

% Computation of the test statistic
n = length(U);
SUM1 = 0; SUM2 = 0;  
for i=1:n
    SUM2 = SUM2 + U(i)^2;
    for j=1:n
        SUM1 = SUM1 + max(U(i),U(j));
    end
end
Sn = SUM2 - (SUM1/n) + (n/3);

% Estimation of the P-value
S = zeros(1,M);
parfor k=1:M
    SUM1 = 0; SUM2 = 0;
    R = unifrnd(0,1,n);
    for i=1:n
        SUM2 = SUM2 + R(i)^2;
        for j=1:n
            SUM1 = SUM1 + max(R(i),R(j));
        end
    end
    S(k) = SUM2 - (SUM1/n) + (n/3);
end
PV = sum(S>Sn) / M;
Test = (PV < 0.05);