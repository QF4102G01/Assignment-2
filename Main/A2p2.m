%Group members: Chen Penghao, Wang Zexin
%Group number: G01

% Initialize values for FDM explicit III
S0 = 9.8;
X = 9;
T = 0.25;
r = 0.001;
q = 0.01;
sigma = 0.15;
N = 2900;
ds = 0.05;

% Execute FDM explicit III
run1 = explicitIIIFDM(S0, X, r, T, q, sigma, N, ds);

% Determine a lower bound for N = T/deltaT
Smax = 4 * X;
I = floor(Smax / ds);
N_lb = T * (sigma ^ 2 * I ^ 2 + r * I);

% Rerun explicit III FDM with the value of lower bound N
run2 = explicitIIIFDM(S0, X, r, T, q, sigma, N_lb, ds);

% Rerun explicit III FDM and lower N gradually and locate the cut-off value
% of N where the estimate loses all its significant figures
N_co = 2916;
run3 = explicitIIIFDM(S0, X, r, T, q, sigma, N_co, ds);

% Run explicit III FDM on American vanilla call option
run4 = explicitIIIFDMAmericanCall(S0, X, r, T, q, sigma, N, ds);