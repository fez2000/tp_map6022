function ScatterPlotRanksNormalised(X,distr, names)

% "Pairwise scatterplots, standardized rank scatterplots & histograms"
% Input  -> X: n x d data matrix
% Output -> d x d matrix of figures structured as follows:
%             * Diagonal: Histograms
%             * Upper triangle: Scatterplots of the original data
%             * Lower triangle: Scatterplots of the standardized ranks
% Necessary procedure: VectorOfRanks

[n,d] = size(X); U = VectorNormalised(VectorOfRanks(X)) / (n+1);

% Diagonal: Histograms
K = 2*(1+3.322*log10(n)); % Number of classes according to Sturge's rule
for i=1:d
    subplot(d,d,(d+1)*i-d);
    disp(round(K))
    histfit(X(:,i), round(abs(K)),distr(i)); 
    title(strcat(["Hist et lois  ",distr(i), "Skew=",skewness(X(:,i))]),'fontsize',6);
    ylabel(names(i),"FontSize",6,"FontName","arial")
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor','k','EdgeColor','w')
end
% Upper triangle: Scatterplots of the original observations
for i=1:d-1
    for j=i+1:d
        POS = (i-1)*d + j;
        subplot(d,d,POS);
        plot(X(:,i),X(:,j),'.k');
       
        ylabel(names(i),"FontSize",6,"FontName","arial")
        xlabel(names(j),"FontSize",6,"FontName","arial")
    end
end
% lower triangle: Scatterplots of the standardized ranks
for i=1:d-1
    for j=i+1:d        
        POS = (j-1)*d + i;
        subplot(d,d,POS);
        plot(U(:,i),U(:,j),'.k'); 
        title(strcat([names(i),' et ', names(j),' normalisé']),'fontsize',6);
    end
end
hold off;