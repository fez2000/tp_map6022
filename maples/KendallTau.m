function Tau_n = KendallTau(X)

% "Sample value of Kendall's measure of association"
% Input  -> X: n x 2 data matrix
% Output -> Tau_n: Empirical Kendall's tau

n = length(X);
SUM = 0;
for j=1:n
    for k=1:n
        SUM = SUM + prod(double(X(j,:)<X(k,:)));
    end
end
Tau_n = 4*SUM/(n*(n-1)) - 1;
