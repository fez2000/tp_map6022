function V = UniformSphere(n,d)

% "Simulated vectors uniformly distributed on the hypersphere"
% Input  -> n: Sample size
%           d: Dimension
% Output -> V: n x d matrix of the generated vectors

mu = zeros(d,1); Sigma = eye(d);

Z = mvnrnd(mu,Sigma,n);
T = sqrt(diag(Z*Z.'));

V = zeros(n,d);
for j=1:d
    V(:,j) = Z(:,j) ./ T;
end