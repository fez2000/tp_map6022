%Loi Normale
%function E = EchantillonageSn()
Sizes = [10,100,500,1000,100000,1000000];
Alphas = [1,2,3,4,5,6,7,8,9,10];

h = @(x)(4*x*(1-x));
ns = 6;
na = 10;
BI = zeros(na,ns);
EQM = zeros(na,ns);
Experices = zeros(na,ns);
for ia = 1:na
    for is=1:ns
        alpha = 2; beta = 1; n = Sizes(is);
        W = betarnd(alpha,beta,n,1);
        Y = zeros(n,1);
        for i=1:n
            disp(h(W(i)));
            Y(i) = h(W(i));
        end
        mu = mean(Y);
        BI(ia,is) = mu - h(0.5);
        EQM(ia,is) = var(Y) + BI(is)^2;
    end
end
%M = 100; n = 10; mu = 0; sigma = 1;
%NmuSn = zeros(M,1);
%NSn = zeros(M,1);
%Sn = zeros(M,1);
%W = normrnd(mu,sigma,n,1);
%sd = std(W)
%p = pdf("Nakagami",W,mean(W),1)
%for i=1:M
%    W = normrnd(mu,sigma,n,1);
%    sd = std(W)
%    NmuSn(i) = pdf("Nakagami",sd,mean(W),1);
%    Sn(i) = std(W);
%    f = pdf("Nakagami",sd,n-1/2,1);
    %hist(f)
%    NSn(i) = f
%end
%figure;
%subplot(2,1,1)
%histogram(Sn);
%title("Sn");
%subplot(2,1,2)
%histogram(NSn)
%title("Desitite Sn");
%end