function ans = Saison (path)
 
 arguments 
    path="datasets\Saison20232024.mat"; 
 end
 x = input("input file path, defauld is: datasets\Saison20232024.mat");
 [Mtilde,M,mu,Sigma_chapeau,R_chapeau, sd] = DefaultStats(x);
 disp(["La moyenne ",mu, "Sigma Chapeau",Sigma_chapeau, "R chapeau ",R_chapeau, "Standart deviation", sd]);
 disp("Work space save at outputs\data\workspaceSaision.mat");
 disp("images save at outputs\figures");
end