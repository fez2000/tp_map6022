function [SK,SP, SPT] = Test_lois_variances_asymptotiques(C, M,d,kendallTauC)
    arguments
        C = 1; % la lois 
        M = 800; % nombre de test
        d = 25; % dimension de l'echantillonage
        kendallTauC = 0; % correlation entre les donnees
    end
   
    
    [SK,SP, SPT] = Kendall_Sperman_asymptotiques(C, M,d,kendallTauC);%generation des donnees
    fig = figure; %ploting surface
    ScatterPlotNormalised([SK, SPT],["Normal","Normal"],["X","Y"]);
   

    savefig("outputs/figures/ScatterPlotAdequations_Test_lois_variances_asymptotiques.fig")
    saveas(fig,"outputs/figures/ScatterPlotAdequations_Test_lois_variances_asymptotiques.png","png")
end
