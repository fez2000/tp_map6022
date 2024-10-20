function [Mtilde,M,mu,Sigma_chapeau,R_chapeau, sd] = Saison()
    load("datasets\Saison20232024.mat");
    
    titles = ["Temps De Jeuen Sec", "Buts Esprs", "Lancers"];
    %calcule de la moyenne 
    M = Saison20232024;
    mu = mean(M);
    
    %calcule de la variance
    
    % Obtention de la matrice des variances-covariances
    Sigma_chapeau = cov(M);
    % Obtention de la matrice des corrélations
    R_chapeau = corr(M);
    sd = std(M);
    % Standardisation
    Mtilde = (M-mu) ./ std(M);
    
    
    
    numCols = 3;
    figure;
    for i=1:numCols
        subplot(2,2,i);
        histogram(Mtilde(:,i));
        title(titles(i));
    end
    savefig("outputs/figures/HistogramsDistributions.fig")
    close;
    
    % Nuage de points par paires 
    
    figure;
    [b,ax] = plotmatrix(Mtilde,Mtilde);
    for i=1:numCols
        ylabel (ax(i,1),titles(i));
        xlabel (ax(numCols,i),titles(i));
    end
    savefig("outputs/figures/ScratterPlot.fig")
    close;
    save outputs\data\workspaceSaision.mat Sigma_chapeau titles mu R_chapeau sd 
    exit(0);
end
