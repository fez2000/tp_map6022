function c = BivCopula_Density(C,TauC,u,v)

% "Copula density of a given copula in List#1b"
% Input  -> C: From List#1b
%           TauC: Kendall's tau
%           u,v: Function's arguments
% Output -> c: Value of c at (u,v)

theta = BivCopula_InversionKendall(C,TauC);

% Normal
if C(1)==1
    if abs(theta)>1
        c = 0.00000001;
    else
        x = norminv(u); y = norminv(v);
        c = mvnpdf([x,y],[0 0],[1 theta;theta 1]) / (normpdf(x)*normpdf(y));
    end

% Student
elseif C(1)==3
    df = C(2);
    if abs(theta)>1
        c = 0.00000001;
    else
        x = tinv(u,df); y = tinv(v,df);
        c = mvtpdf([x,y],[1 theta;theta 1],df) / (tpdf(x,df)*tpdf(y,df));
    end

% Laplace
elseif C(1)==5
    beta = C(2);
    Fu = @(r)(Laplace_cdf(r,beta)-u);
    Fv = @(r)(Laplace_cdf(r,beta)-v); 
    w = [fzero(Fu,0) fzero(Fv,0)];
    
    sigma = [1 theta; theta 1];
    L1 = sqrt(w/sigma*w.'); L2 = 1/sqrt(det(sigma));
    
    d = 2; a = beta - d/2;
    A = 2*2^(-a/2)*L2*L1^a*besselk(a,sqrt(2)*L1)/((2*pi)^(d/2)*gamma(beta));
    
    a = beta - 1/2; 
    B1 = 2*2^(-a/2)*(abs(w(1)))^a*besselk(a,sqrt(2)*abs(w(1)))/(sqrt(2*pi)*gamma(beta));
    B2 = 2*2^(-a/2)*(abs(w(2)))^a*besselk(a,sqrt(2)*abs(w(2)))/(sqrt(2*pi)*gamma(beta));
    c = A / (B1*B2);
      
% Plackett
elseif C(1)==9
    if abs(theta)>1
        c = 0.00000001;
    elseif theta==0
        c = 1;
    else
        T = (1+theta)/(1-theta);
        L = 1 + (T-1)*(u+v);
        A = T*(L + 2*(1-T)*u*v);
        B = (u+v-1)^2 + 2*T*(u*(1-u)+v*(1-v)) + T^2*(u-v)^2;
        c = A / B^(3/2); 
    end
  
% Clayton 
elseif C(1)==11
    if theta<0
        c = 0.00000001;
    elseif (theta==0)
        c = 1;
    else
        c = (theta+1)*(u*v)^(-theta-1)*( u^(-theta) + v^(-theta) - 1 )^(-(1/theta)-2);
    end
 
% Frank
elseif C(1)==13
    if theta<0
        c = 0.00000001;
    elseif theta==0
        c = 1;
    else
        r = exp(-theta); s = exp(-theta*u); t = exp(-theta*v);
        c = theta*(1-r)*s*t/(r-s-t+s*t)^2;
    end
    
% Gumbel
elseif C(1)==15
    if (theta<0) || (theta>=1)
        c = 0.00000001;
    else
        % A finite-difference is used to increase numerical stability
        epsilon = 0.00001;
        c1 = BivCopula_Copula(15,theta,u+epsilon,v+epsilon) + BivCopula_Copula(15,theta,u,v);
        c2 = BivCopula_Copula(15,theta,u+epsilon,v) + BivCopula_Copula(15,theta,u,v+epsilon);
        c = (c1-c2)/epsilon^2;
    end
    
% Farlie-Gumbel-Morgenstern
elseif C(1)==17
    if abs(theta)>1
        c = 0.00000001;
    else
        c = 1 + theta*(1-2*u)*(1-2*v);
    end

% Squared copulas of the Normal (centred Chi-square), Student (Fisher), 
% Laplace, Pearson II, Plackett, Clayton, Frank, Gumbel
elseif C(1)==2 || C(1)==4 || C(1)==6 || C(1)==8 || C(1)==10 || C(1)==12 || C(1)==14 || C(1)==16
    C(1) = C(1)-1;
    u1 = (1+u)/2; u2 = (1-u)/2; v1 = (1+v)/2; v2 = (1-v)/2;
    c1 = BivCopula_Density(C,theta,u1,v1); 
    c2 = BivCopula_Density(C,theta,u1,v2);
    c3 = BivCopula_Density(C,theta,u2,v2); 
    c4 = BivCopula_Density(C,theta,u2,v1);
    c = (c1+c2+c3+c4)/4;     
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F = Laplace_cdf(x,beta)

f = @(r)(normcdf(x./sqrt(r)).*gampdf(r,beta,1));
F = integral(f,0,Inf);
end