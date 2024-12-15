function R = VectorNormalised(X) 

% "Vectors of ranks from multivariate data"
% Input  -> X: n x d data matrix
% Output -> R: n x d matrix of ranks computed for each column

[n,d] = size(X);
R = zeros(n,d);
for j=1:d
    R(:,j) = normalize(X(:,j));
end