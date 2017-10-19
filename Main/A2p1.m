rhos = [1 1/2 1/5];
Ns = [50 100 200 400];

for rho = rhos
    for N = Ns
        disp(['For rho = ', num2str(rho), ' and N = ', num2str(N), '. The option price is : ']);
        FSGMAmericanFixedStrikeAsianPut(0.25, 0.25, 100, 0.4, 0.01, 95, 0.1, 100, N, rho)
    end
end

Ns = 50 : 50 : 500;
for N = Ns
    FSGMAmericanFixedStrikeLookbackPut(0.25, 0.25, 1, 0.4, 0.01, 0.97, 0.1, 0.95, N)
end
for N = Ns
    FSGMAmericanFixedStrikeLookbackPut(0.25, 0.25, 1, 0.4, 0.01, 0.57,  0.1, 0.95, N)
end