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
disp(['Executing explicit scheme III with N = ', num2str(N)]);
run1 = explicitIIIFDM(S0, X, r, T, q, sigma, N, ds);
disp(['The european vanilla call option price for this run is ', num2str(run1)]);

% Compare against the exact value
bench = EuropeanVanillaCall(S0, q, X, T, r, sigma);
disp(['The value of the call option using the formulae is ', num2str(bench)]);

% Determine a lower bound for N = T/deltaT
Smax = 4 * X;
I = floor(Smax / ds);
N_lb = round(T * (sigma ^ 2 * I ^ 2 + r * I));

% Rerun explicit III FDM with the value of lower bound N
disp(['Executing explicit scheme III with N = ', num2str(N_lb)]);
run2 = explicitIIIFDM(S0, X, r, T, q, sigma, N_lb, ds);
disp(['The european vanilla call option price for this run is ', num2str(run2)]);

% Rerun explicit III FDM and lower N gradually and locate the cut-off value
% of N where the estimate loses all its significant figures
disp(['Executing explicit scheme III with N = ', num2str(N_co)]);
N_co = 2900;
run3 = explicitIIIFDM(S0, X, r, T, q, sigma, N_co, ds);
disp(['The european vanilla call option price for this run is ', num2str(run3)]);

% Run explicit III FDM on American vanilla call option
disp(['Executing explicit scheme III on American vanilla call with N = ', num2str(N)]);
run4 = explicitIIIFDMAmericanCall(S0, X, r, T, q, sigma, N, ds);
disp(['The american vanilla call option price for this run is ', num2str(run4)]);