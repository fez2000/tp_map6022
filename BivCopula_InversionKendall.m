function theta = BivCopula_InversionKendall(C,TauC)

% "Parameter of a copula in List#1a corresponding to a given Kendall's tau"
% Input  -> C: From List#1a
%           TauC: Kendall's tau
% Output -> theta: Copula parameter
% Necessary tables: Kendall_Plackett.txt, Kendall_Frank.txt, Kendall_HuslerReiss.txt

% Normal, Student, Laplace & Pearson Type II
if C(1)==1 || C(1)==3 || C(1)==5 || C(1)==7
    theta = sin((pi/2)*TauC);

% Centred Chi-square
elseif C(1)==2
    TauC = max(0,TauC);
    theta = sin((pi/2)*sqrt(TauC));
    
% Fisher, Squared-Laplace, Squared-Pearson type II
elseif C(1)==4 || C(1)==6 || C(1)==8
    rng('shuffle');
    S = 250;
    f = @(x)(Kendall_SquaredElliptical(x,C,S)-TauC)^2;
    LB = .01; UB = .99; x0 = (LB+UB)/2;
    options = optimset('maxiter',50);
    theta = fminsearchbnd(f,x0,LB,UB,options);  
    
% Plackett
elseif C(1)==9
    f = @(x)(abs(Kendall_Plackett(x)-TauC));
    options = optimset('maxiter',25);
    LB = -.99; UB = .99; x0 = (LB+UB)/2;
    theta = fminsearchbnd(f,x0,LB,UB,options);
    
% Clayton
elseif C(1)==11
    theta = 2*TauC/(1-TauC);
    
% Frank
elseif C(1)==13
    f = @(x)(Kendall_Frank(x)-TauC);
    theta = fzero(f,.5);
    
% Gumbel
elseif C(1)==15
    theta = max(0,TauC);

% Squared-Plackett, Squared-Clayton, Squared-Frank, Squared-Gumbel
elseif C(1)==10 || C(1)==12 || C(1)==14 || C(1)==16   
    rng('shuffle'); S = 250; D = C(1)-1;
    f = @(x)(Kendall_SquaredArch(x,D,S)-TauC)^2;
    LB = .01; UB = .99; x0 = (LB+UB)/2;
    options = optimset('maxiter',50);
    TauTilde = fminsearchbnd(f,x0,LB,UB,options);
    theta = BivCopula_InversionKendall(D,TauTilde);

% Farlie-Gumbel-Morgenstern
elseif C(1)==17
    if TauC < -2/9
        theta = -1;
    elseif TauC > 2/9
        theta = 1;
    else
        theta = 9*TauC/2;
    end 
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Tau = Kendall_SquaredElliptical(x,C,S)

if C(1)==4
    Y = mvtrnd([1 x;x 1],C(2),S);
    Tau = KendallTau(Y.^2);
elseif C(1)==6
    Z = mvnrnd([0 0],[1 x;x 1],S); 
    R = gamrnd(C(2),1,S,1);
    Y = [ diag(sqrt(R)*Z(:,1).'), diag(sqrt(R)*Z(:,2).') ];
    Tau = KendallTau(Y.^2);
elseif C(1)==8
    A = chol([1 x;x 1]).';
    W = zeros(S,2);
    for i=1:S
        G = sqrt(betarnd(1,C(2)+1,1));
        V = UniformSphere(1,2);
        W(i,:) = G*A*V.';
    end
    U = (1+sign(W).*betacdf(W.^2,1/2,C(2)+3/2)) / 2;
    Tau = KendallTau(abs(2*U-1));
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Tau = Kendall_Plackett(theta)

alpha = 2*theta/(1-theta);

A = @(v)((alpha+1).*v - alpha.*v.^2);
B = @(v)(integral(@(s)((alpha*(alpha+1).*(1-v).^2.*s.^2 + (alpha+1).*v - alpha.*v^2 + (alpha+1).*(1-v).*s).^(-1)),0,1));

T = 200;
K = zeros(T,1);
for j=1:T
    v = (j-.5) / T;
    K(j) = v + (1-v).*A(v).*B(v);
end
Tau = 3 - 4*mean(K);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Tau = Kendall_Frank(theta)

f = @(t)(t./(exp(t)-1));
Tau = 1 - (4/theta) + (4/theta^2)*integral(f,0,theta);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Tau = Kendall_SquaredArch(x,D,S)

U = BivCopula_Simulation(S,D,x);
Y = abs(2*U-1);
Tau = KendallTau(Y);
end