%Group members: Chen Penghao, Wang Zexin
%Group number: G01

% Initialize values for FDM explicit III
S0 = 9.8;
X = 9;
T = 0.25;
r = 0.001;
q = 0.01;
sigma = 0.15;
dt = 0.01;
ds = 0.05;
N = T/dt;

% Execute FDM explicit III
disp(['Executing explicit scheme III with N = ', num2str(N)]);
run1 = explicitIIIFDM(S0, X, r, T, q, sigma, N, ds, false, true);
disp(['The approximated european vanilla call option price for this run is ', num2str(run1)]);

% Determine a lower bound for N = T/deltaT
Smax = 4 * X;
I = floor(Smax / ds);
N_lb = ceil(T * (sigma ^ 2 * I ^ 2 + r * I));

% Rerun explicit III FDM with the value of lower bound N
disp(['Executing explicit scheme III with N = ', num2str(N_lb)]);
run2 = explicitIIIFDM(S0, X, r, T, q, sigma, N_lb, ds, true, false);
disp(['The approximated european vanilla call option price for this run is ', num2str(run2)]);

% Rerun explicit III FDM and lower N gradually and locate the cut-off value
% of N where the estimate loses all its significant figures

check = run2;
for n = N_lb : -1 : 2440
    run3 = explicitIIIFDM(S0, X, r, T, q, sigma, n, ds, false, false);
    if round(run3,5) ~= round(check, 5)
        disp(['Executing explicit scheme III with N = ', num2str(n)]);
        disp(['The approximated european vanilla call option price for this run is ', num2str(run3)]);
    end
end

N_cutoff = 2453;

% Run explicit III FDM on American vanilla call option
disp(['Executing explicit scheme III on American vanilla call with N = ', num2str(N)]);
run4 = explicitIIIFDMAmericanCall(S0, X, r, T, q, sigma, N_lb, ds);
disp(['The approximated american vanilla call option price for this run is ', num2str(run4)]);