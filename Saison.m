function ans = Saison (path)
 
 arguments 
    path="datasets\Saison20232024.mat"; 
 end
 x = input("input file path, defauld is: datasets\Saison20232024.mat");
 if isempty(x)
    x = path;
 end
 [Mtilde,M,mu,Sigma_chapeau,R_chapeau, sd] = DefaultStats(x);
 out = fprintf('La moyenne %d, Sigma Chapeau %d, R chapeau %d, Standart deviation %d\n',mu, Sigma_chapeau, R_chapeau, sd);
 disp(out);
 disp("Work space save at outputs\data\workspaceSaision.mat");
 disp("images save at outputs\figures");
end