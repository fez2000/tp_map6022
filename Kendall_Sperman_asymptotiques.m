function [SK,SP, SPT] = Kendall_Sperman_asymptotiques(C, M,d,kendallTauC)


    SK = zeros(M,1);
    SPT = zeros(M,1);
    SP= zeros(M,1);
% Générer M échantillons de taille n à partir d'une distribution normale standard. 

    for n=1:M
        U = BivCopula_Simulation(d,C,kendallTauC);
        SK(n) = sqrt(n) * KendallTau(U) *3/2;
        s = SpearmanRho(U);
        SP(n) = s; 
        SPT(n) = sqrt(n) * s;
    end
    mu = mean(SP);
    for n=1:M
       
        SPT(n) = sqrt(n) * (SP(n)-mu);
    end

end
