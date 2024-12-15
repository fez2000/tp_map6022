function U = BivCopula_Simulation(n,C,TauC)

% "Simulated pairs from a given copula in List#1a"
% Input  -> n: Sample size
%           C: From List#1a
%           TauC: Kendall's tau
% Output -> U: n x 2 matrix of generated pairs
% Necessary procedures: BivCopula_InversionKendall

if TauC==0
    U = rand(n,2);
elseif TauC==1
    X = rand(n,1);
    U = repmat(X,1,2);
elseif TauC==-1
    X = rand(n,1);
    U(:,1) = X; U(:,2) = 1-X;
else
    U = zeros(n,2);
    theta = BivCopula_InversionKendall(C,TauC);

    % Normal & Centred Chi-square
    if C(1)==1 || C(1)==2
        mu = [0 0]; Sigma = [1 theta; theta 1];
        X = mvnrnd(mu,Sigma,n);
        U = normcdf(X);
        if C(1)==2
            U = abs(2*U-1);
        end

    % Student & Fisher
    elseif C(1)==3 || C(1)==4
        X = mvtrnd([1 theta;theta 1],C(2),n);
        U = tcdf(X,C(2));
        if C(1)==4
            U = abs(2*U-1);
        end
        
    % Laplace & Squared-Laplace
    elseif C(1)==5 || C(1)==6
        Z = mvnrnd([0 0],[1 theta;theta 1],n); 
        R = gamrnd(C(2),1,n,1);
        Y = sqrt(R).*Z;
        for i=1:n
            U(i,1) = Laplace_cdf(Y(i,1),C(2));
            U(i,2) = Laplace_cdf(Y(i,2),C(2));
        end
        if C(1)==6
            U = abs(2*U-1);
        end
    
    % Pearson type II & Squared-Pearson type II
    elseif C(1)==7 || C(1)==8
        A = chol([1 theta; theta 1]).'; W = zeros(n,2);
        for i=1:n
            G = sqrt(betarnd(1,C(2)+1,1));
            V = UniformSphere(1,2);
            W(i,:) = G*A*V.';
        end
        U = (1+sign(W).*betacdf(W.^2,1/2,C(2)+3/2)) / 2;
        if C(1)==8
            U = abs(2*U-1);
        end
    
    % Plackett & Squared-Plackett
    elseif C(1)==9 || C(1)==10
        for i=1:n
            V = rand; T = rand;
            a = T*(1-T); 
            b = ( 1 - theta^2 + 4*a*theta^2 ) / (1-theta)^2;
            c = ( 2*a*( (1+theta)^2*V + (1-theta)^2*(1-V) ) + (1-theta^2)*(1-2*a) ) / (1-theta)^2;
            d = sqrt( (1+theta)*( 1 - theta^2 + 16*a*V*(1-V)*theta^2 ) ) / (1-theta)^(3/2);
            U(i,1) = V; U(i,2) = (c-(1-2*T)*d) / (2*b);
        end
        if C(1)==10
            U = abs(2*U-1);
        end

    % Clayton & Squared-Clayton
    elseif C(1)==11 || C(1)==12
        for i=1:n
            V = rand; T = rand;
            a = V^(-theta);
            b = T^(-theta/(theta+1)) - 1;
            U(i,1) = V; U(i,2) = (a*b+1)^(-1/theta);              
        end
        if C(1)==12
            U = abs(2*U-1);
        end
       
    % Frank & Squared-Frank
    elseif C(1)==13 || C(1)==14
        for i=1:n
            V = rand; T = rand;
            a = T*exp(-theta) + (1-T)*exp(-theta*V);
            b = T + (1-T)*exp(-theta*V);
            U(i,1) = V; U(i,2) = log(b/a) / theta;
        end
        if C(1)==14
            U = abs(2*U-1);
        end
       
    % Gumbel & Squared-Gumbel
    elseif C(1)==15 || C(1)==16
        for i=1:n
            T0 = rand; T1 = rand; T2 = rand;
            Z = (T0^(1-theta))/(T0^(1-theta)+(1.0-T0)^(1-theta));
            if (rand < 1-theta) 
                W = T1*T2;
            else
                W = T2;
            end
            R = ((1-Z)^(1/(1-theta)) + Z^(1/(1-theta)))^(1-theta);
            U(i,1) = W^(Z/R); U(i,2) = W^((1.0-Z)/R);
        end
        if C(1)==16
            U = abs(2*U-1);
        end
    
    % Farlie-Gumbel-Morgenstern
    elseif C(1)==17
        if theta == 1
            U(:,1) = rand(n,1); U(:,2) = U(:,1);
        else
            for i=1:n
                V = rand; T = rand;
                b = theta*(1-2*V);
                R = sqrt((1+b)^2 - 4*T*b);
            U(i,1) = V; U(i,2) = 2*T/(1+b+R);
            end
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F = Laplace_cdf(x,beta)

f = @(r)(normcdf(x./sqrt(r)).*gampdf(r,beta,1));
F = integral(f,0,Inf);
end