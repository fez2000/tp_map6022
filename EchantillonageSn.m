%Loi Normale
%function E = EchantillonageSn()
alpha = 2; beta = 1; n = 1000;
W = betarnd(alpha,beta,n,1);
h = @(x)(4*x(1-x)) 
for i=1:n
    h(W(i))
%for i=1:M
%    W = normrnd(mu,sigma,n,1);
%    sd = std(W)
%    NmuSn(i) = pdf("Nakagami",sd,mean(W),1);
%    Sn(i) = std(W);
%    f = pdf("Nakagami",sd,n-1/2,1);
    %hist(f)
%    NSn(i) = f
%end
figure;
subplot(2,1,1)
histogram(Sn);
title("Sn");
subplot(2,1,2)
histogram(NSn)
title("Desitite Sn");
%end