function ans = Variances_asymptotiques(C, M,d,kendallTauC)
    arguments
        C = 15;
        M = 800;
        d = 25;
        kendallTauC = 0;
    end
   

    [SK,SP, SPT] = Kendall_Sperman_asymptotiques(C, M,d,kendallTauC);
  
    ans = [var(SK),var(SPT)] ;

end