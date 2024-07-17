function v = bigstressfunc(t,Y,f,d)

% defining parameters
R_0 = 36;
E_GO = 0.06;
% S_I = 0.01; %0.03;
sig = 1.8;
alp = 20000;
K = 18;
d_0 = 0.0025;
r_1 = 0.000035;
r_2 = 0.0000001;
G_e = 5.83;
rho = 0;

kcd=1;
kad=10;
kcr=0.05;
krd=0.3;
koi1=0.2;
koi2=0.2;
k_hpa = 0.01;
kc = 2.2;

%parameters for updated model

% Default Parameter Values
% S_0 = 0.01;              % ml/(µU·day), example value
% V = 0.01;               % l
V = 0.01 * 500;
k = 700;                % day^-1
tau_gamma = 2.14;       % days
tau_sigma = 10.71;      % days
tau_beta = 42.85;       % days
% tau_IS = 0.1;           % days, example value

% Auxiliary Function Parameters
kM = 2;                 % unitless
alpha_M = 150;          % mg/dl
kISR = 2;               % unitless
alpha_ISR = 1.2;        % unitless
P_max = 4.55;           % day^-1
kP = 4;                 % unitless
alpha_P = 41.77;        % unitless
A_max = 3.11;           % day^-1
kA = 6;                 % unitless
alpha_A = 0.44;         % unitless
A_b = 0.8;              % day^-1
gamma_max = 0.2;        % unitless
gamma_s = 99.9;         % unitless
gamma_n = 1;            % unitless
gamma_0 = 0.1;          % unitless
sigma_ISRmax = 867.6;   % µU·m/mg·dl
sigma_ISRs = 0.1;       % unitless
sigma_ISRn = 0.1;       % unitless
sigma_ISRk = 1;         % unitless
sigma_b = 3;            % µU·m/mg·dl
sigma_Mmax = 1;         % unitless
sigma_Ms = 0.2;         % unitless
sigma_Mn = 0.02;        % unitless
sigma_Mk = 0.2;         % unitless

S_0 = 0.022/24;
tau_S_I = 5.1*24;

function M_val = M(G)
    M_val = ((G^kM) / (alpha_M^kM + G^kM));
end

function ISR_val = ISR(G, gamma, sigma)
    ISR_val = sigma * (((M(G)+ gamma)^kISR) / (alpha_ISR^kISR + ((M(G) + gamma)^kISR)));
end

function P_val = P(ISR)
    P_val = (P_max) * ((ISR^kP) / (alpha_P^kP + ISR^kP));
end

function A_val = A(M)
    A_val = A_max * ((M^kA) / (alpha_A^kA + M^kA));
end

function gamma_inf_val = gamma_inf(G)
    gamma_inf_val = ((gamma_max) / 1 + exp((G - gamma_s) / gamma_n)) - gamma_0;
end

function sigma_inf_val = sigma_inf(ISR, M, G, gamma, sigma)
    sigma_inf_val = sigma_ISR_inf(ISR, G, gamma, sigma) * sigma_Minf(M, G) + sigma_b;
end

function sigma_ISR_inf_val = sigma_ISR_inf(ISR, G, gamma, sigma)
    sigma_ISR_inf_val = (sigma_ISRmax) / (1 + sigma_ISRk * exp(-(ISR_sigma(M_sigma(G), gamma, sigma) - sigma_ISRs) / sigma_ISRs));
end

function sigma_Minf_val = sigma_Minf(M, G)
    sigma_Minf_val = 1 - ((sigma_Mmax) / 1 + sigma_Mk * exp(-(M_sigma(G) - sigma_Ms) / sigma_Ms));
end

function ISR_sigma_val = ISR_sigma(M_sigma, gamma, sigma)
    ISR_sigma_val = ISR(M_sigma, gamma, sigma);
end

function M_sigma_val = M_sigma(G)
    G_sigmas = 0;
    M_sigma_val = M(G - G_sigmas);
end

G = Y(1);
I = Y(2);
B = Y(3);
c = Y(4);
a = Y(5);
r = Y(6);
o = Y(7);  
gamma = Y(8);
sigma = Y(9);
S_I = Y(10);


v(1) =  (kc*o*G)+ R_0  - (E_GO+S_I*I)*G;
v(2) = (1-rho)*((B*(sig)*G^2)/(alp +G^2))-(K)*I; %from Topp
% v(2) = (B / V) * ISR(M(G), gamma, sigma) - k*I;
% v(3) = ((-d_0 + r_1*G - r_2*G^2) * B); %from Topp
% v(3) = ((P(ISR(G, gamma, sigma)) - A(M(G)))*B) / (tau_beta); % new beta e
v(3) =  (((P(ISR(G, gamma, sigma)) - A(M(G)))*B) / (tau_beta)) / 200;
v(4) = (1+f)/(1+o/koi1) - kcd*c;
v(5) = c/(1+(o+d)*r/koi2) - kad*a;
% v(5) = c/(1+((o*r)/koi2)) - kad*a;
v(6) = ((o+d)*r)*((o+d)*r)/(k_hpa+((o+d)*r)*((o+d)*r)) + kcr - krd*r;
% v(6) = ((o*r)^2)/(k + (o*r)^2) + kcr + krd*r;
v(7) = a - o;

v(4) = (1+f)/(1+o/koi1) - kcd*c;
v(5) = c/(1+(o+d)*r/koi2) - kad*a;
v(6) =((o+d)*r)*((o+d)*r)/(k+((o+d)*r)*((o+d)*r)) + kcr - krd*r;
v(7) = a - o;

% v(8) = (gamma_inf(ISR(G, gamma, sigma)) - gamma)/tau_gamma; % dgamma/dt
v(8) = (gamma_inf(G) - gamma)/tau_gamma;
% v(9) = (sigma_inf(ISR(G, gamma, sigma), M(G), G, gamma, sigma) - sigma)/tau_sigma; % dsigma/dt
v(9) = ((sigma_inf(ISR(G, gamma, sigma), M(G), G, gamma, sigma) - sigma)/tau_sigma) / 23.3;
v(10) = (S_0 - S_I) / tau_S_I;
v = v';
end