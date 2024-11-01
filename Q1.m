rng default;  % For reproducibility
function S = QNakagami(M,n)


    S = zeros(M,1);

% Générer M échantillons de taille n à partir d'une distribution normale standard. 

    for i=1:M
        X = normrnd(0,1,n,1);
        S(i) = std (X);
    end
end
% Distribution de Nakagami aux écarts types de l'échantillon. 
S = QNakagami(1000, 10);
fig = figure; %ploting surface
histfit(S,100,'nakagami');
savefig("outputs/figures/HistogramsDistributions_normal_nakagami.fig")
saveas(fig,"outputs/figures/HistogramsDistributions_normal_nakagami.png","png")
close;
%la loi de Nakagami est moins robuste lorsque les données proviennent d'une distribution proche de la normale

%%%%%%%% b) Robustesse sur la loi de student
%%%%%%Pour nu=[1:3,6,10]
function S = QStudent(M,n)
     S = zeros(M,1);
    for i=1:M
        W = trnd(5,n,1);
        S(i) = std (W);
    end
end
S = QStudent(1000, 10);
fig = figure;
histfit(S,100,'nakagami');
savefig("outputs/figures/HistogramsDistributions_student_nakagami.fig")
saveas(fig,"outputs/figures/HistogramsDistributions_student_nakagami.png","png")
close;
% Apres plusieurs reprises avec les valeurs de nu on constate que plus nu est grand, plus la loi de Student se rapproche d'une loi normale.
% L ecart entre la courbe de desnsit2 et l histogramme se remarque moins a chaque valeur de nu.



%%%%%%%% c) Robustesse sur la loi gamma
%Pour alpha=[1:5]
%Loi Gamma (alpha, 1)
function S = QGamma(M,n)
    S = zeros(M,1);
    alpha = 5;
    beta = 2;
    for i=1:M
        W = gamrnd(alpha, beta,n,1);
        S(i) = std (W);
    end
end

%Plus alpha est grand, plus la loi gamma se rapproche d'une loi normale et est robuste.
fig = figure;
S = QGamma(1000, 10);
histfit(S,100,'nakagami');
savefig("outputs/figures/HistogramsDistributions_gamma_nakagami.fig");
saveas(fig,"outputs/figures/HistogramsDistributions_gamma_nakagami.png","png");
close;