function ScatterPlotRanks(X)

% "Pairwise scatterplots, standardized rank scatterplots & histograms"
% Input  -> X: n x d data matrix
% Output -> d x d matrix of figures structured as follows:
%             * Diagonal: Histograms
%             * Upper triangle: Scatterplots of the original data
%             * Lower triangle: Scatterplots of the standardized ranks
% Necessary procedure: VectorOfRanks

[n,d] = size(X); U = VectorOfRanks(X) / (n+1);

% Diagonal: Histograms
K = 2*(1+3.322*log10(n)); % Number of classes according to Sturge's rule
for i=1:d
    subplot(d,d,(d+1)*i-d);
    hist(X(:,i),K); 
    
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','k','EdgeColor','w')
end
% Upper triangle: Scatterplots of the original observations
for i=1:d-1
    for j=i+1:d
        POS = (i-1)*d + j;
        subplot(d,d,POS);
        plot(X(:,i),X(:,j),'.k');
    end
end
% lower triangle: Scatterplots of the standardized ranks
for i=1:d-1
    for j=i+1:d        
        POS = (j-1)*d + i;
        subplot(d,d,POS);
        plot(U(:,i),U(:,j),'.k');    
    end
end
hold off;