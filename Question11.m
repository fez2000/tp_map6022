%Loi Normale
function [Mtilde,M,mu,Sigma_chapeau,R_chapeau, sd] = EchantillonageSn()
Sizes = [10,100,500,1000,50000,100000,500000,1000000];
Alphas = [1,2,3,4,5,6,7,8,9,10];

h = @(x)(4*x*(1-x));
ns = size(Sizes);
ns = ns(2);
na = size(Alphas);
na = na(2);
BI = zeros(na,ns);
EQM = zeros(na,ns);
E = @(x)(mean(x));
ESB = @(x)(h(0.5));
BE = @(x, reelMu)( E(x) - reelMu );
for ia = 1:na
    for is=1:ns
        alpha = Alphas(ia); beta = 1; n = Sizes(is);
        X = betarnd(alpha,beta,n,1);
        Y = zeros(n,1);
        for i=1:n
            Y(i) = h(X(i));
        end

        BI(ia,is) = BE(Y, ESB(Y));
        EQM(ia,is) = var(Y) + BI(is)^2;
    end
end
end