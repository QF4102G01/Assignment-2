%Group members: Chen Penghao, Wang Zexin
%Group number: G01

% Initialize possible rho values
rhos = [1 1/2 1/5];

% Number of time periods in the lattice
Ns = [50 100 200 400];

% Calculate FSGM American Fixed Strike Asian Put value 
% with the given values of rho and N
for rho = rhos
    for N = Ns
        disp(['For rho = ', num2str(rho), ' and N = ', num2str(N), '. The option price is : ']);
        FSGMAmericanFixedStrikeAsianPut(0.25, 0.25, 100, 0.4, 0.01, 95, 0.1, 100, N, rho)
    end
end

% Calculate FSGM American Fixed Strike Asian Lookback Put value
% with N values ranging from 50 to 500 in increments of 50
Ns = 50 : 50 : 500;
for N = Ns
    FSGMAmericanFixedStrikeLookbackPut(0.25, 0.25, 1, 0.4, 0.01, 0.97, 0.1, 0.95, N)
end

% Calculate FSGM American Fixed Strike Asian Lookback Put value
% by changing running minimum to $0.57
for N = Ns
    FSGMAmericanFixedStrikeLookbackPut(0.25, 0.25, 1, 0.4, 0.01, 0.57, 0.1, 0.95, N)
end