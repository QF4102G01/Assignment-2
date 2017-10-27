%Group members: Chen Penghao, Wang Zexin
%Group number: G01

function [optionValue] = FSGMAmericanFixedStrikeAsianPut(runningTime, tau, S0, sigma, q, runningAverage, r, K, N, rho)
    % Initializing the parameters
    deltaT = tau / N;
    deltaX = sigma * sqrt(deltaT);
    deltaY = rho * deltaX;
    u = exp(deltaX);
    d = exp(-deltaX);
    p = (exp(r * deltaT) - d) / (u-d);
    m = 1 / rho;
    
    % Include the price at t = 0.25 in the Running Average
    runningN = runningTime / deltaT + 1;
    runningAverage = (runningAverage * (runningN - 1) + S0) / runningN;
    
    % Initialize shifts in j and k values
    jshift = 1;
    kshift = N * m + 1;
    
    % Initialize terminal payoff
    v = zeros(N+1, 2 * N * m + 1);
    for j = 0 : N
        for k = -N * m : N * m
            A = runningAverage * exp(k * deltaY);
            v(j+jshift, k+kshift) = max(K - A, 0);
        end
    end
    
    % Backward iterations
    for n = N-1 : -1 : 0
        newV = zeros(N+1, 2 * N * m + 1);
        for j = 0 : n
            for k = -n * m : n * m
                S = S0 * exp((2 * j - n) * deltaX);
                A = runningAverage * exp(k * deltaY);
                
                % Running average for the up and down case
                Au = (A * (n + runningN) + S * u) / (n + runningN + 1);
                Ad = (A * (n + runningN) + S * d) / (n + runningN + 1);
                
                % Obtain the index in the FSGM
                ku = min(floor(log(Au / runningAverage) / deltaY), N * m-1);
                kd = min(floor(log(Ad / runningAverage) / deltaY), N * m-1);
                
                ju = j + 1;
                jd = j;
                
                % Interpolate to get the payoff for up case
                Auup = runningAverage * exp((ku+1) * deltaY);
                Audown = runningAverage * exp(ku * deltaY);
                alphau = (Auup - Au) / (Auup - Audown);
                Vu = (1 - alphau) * v(ju+jshift, ku+1+kshift) + alphau * v(ju+jshift, ku+kshift);
                
                % Interpolate to get the payoff for the down case
                Adup = runningAverage * exp((kd+1) * deltaY);
                Addown = runningAverage * exp(kd * deltaY);
                alphad = (Adup - Ad) / (Adup - Addown);
                Vd = (1 - alphad) * v(jd+jshift, kd+1+kshift) + alphad * v(jd+jshift, kd+kshift);
                
                % Obtain the option value for that point
                newV(j+jshift, k+kshift) = max(exp(-r * deltaT) * (p * Vu + (1-p) * Vd), K - A);
            end
        end
        v = newV;
    end
    
    % Final result
    optionValue = max(v(jshift, kshift), K - runningAverage);
    
end