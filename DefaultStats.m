function [Mtilde,M,mu,Sigma_chapeau,R_chapeau, sd] = DefaultStats(path, titles)
   
    arguments 
        path="datasets\Saison20232024.mat"; 
        titles = ["Temps De Jeuen Sec", "Buts Esprs", "Lancers"];
    end 
    if isfile(path)
        load(path);
    else
    disp([x, "Not Found"]);
    end    
    
    %calcule de la moyenne 
    M = Data;
    mu = mean(M);
    
    %calcule de la variance
    
    % Obtention de la matrice des variances-covariances
    Sigma_chapeau = cov(M);
    % Obtention de la matrice des corrélations
    R_chapeau = corr(M);
    sd = std(M);
    % Standardisation
    Mtilde = (M-mu) ./ std(M);

    numCols = size(M);
    numCols = numCols(2);
    fig =figure;
    for i=1:numCols
        subplot(2,2,i);
        histogram(Mtilde(:,i));
        title(titles(i));
    end
    
    savefig("outputs/figures/HistogramsDistributions.fig")
    saveas(fig,"outputs/figures/HistogramsDistributions.png","png")
    close;
    
    % Nuage de points par paires 
    
    fig = figure;
    [b,ax] = plotmatrix(Mtilde,Mtilde);
    for i=1:numCols
        ylabel (ax(i,1),titles(i));
        xlabel (ax(numCols,i),titles(i));
    end
    savefig("outputs/figures/ScratterPlot.fig")
    saveas(fig,"outputs/figures/ScratterPlot.png","png")
    close;
    save outputs\data\workspaceSaision.mat Sigma_chapeau titles mu R_chapeau sd 
   
end
