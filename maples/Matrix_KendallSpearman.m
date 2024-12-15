function A = Matrix_KendallSpearman(X)

% Matrix of pair-by-paire Kendall and Spearman measures of association"
% Input  -> X: n x d data matrix
% Output -> d x d matrix structured as follows:
%             * Upper triangle: Pairwise Kendall's tau
%             * Lower triangle: Pairwise Spearman's rho
% Necessary procedures: KendallTau, SpearmanRho

d = size(X,2); A = ones(d,d);
for i=1:d-1
    for j=i+1:d
        Y = [X(:,i) X(:,j)];
        A(i,j) = KendallTau(Y);
        A(j,i) = SpearmanRho(Y);
    end
end
    