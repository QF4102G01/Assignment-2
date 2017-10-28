%Group members: Chen Penghao, Wang Zexin
%Group number: G01

function [optionValue] = FSGMAmericanFixedStrikeLookbackPutShoot(runningTime, tau, S0, sigma, q, runningMin, r, K, N)
    % Initializing the parameters
    deltaT = tau / N;
    deltaX = sigma * sqrt(deltaT);
    u = exp(deltaX);
    d = exp(-deltaX);
    p = (exp((r-q) * deltaT) - d) / (u-d);
    
    % Initialize shifts in j and k values
    jshift = 1;
    kshift = N + 1;
    
    % Initialize running minimum for lookback put
    runningMin = min(runningMin, S0);
    
    % Initialize terminal payoff
    % For lookback put, rho is chosen to be 1
    v = zeros(N+1, 2 * N + 1);
    for j = 0 : N
        for k = -N : N
            A = runningMin * exp(k * deltaX);
            v(j+jshift, k+kshift) = max(K - A, 0);
        end
    end
    
    % Backward iterations
    for n = N-1 : -1 : 0
        newV = zeros(N+1, 2 * N + 1);
        for j = 0 : n
            for k = -n : n
                S = S0 * exp((2 * j - n) * deltaX);
                A = runningMin * exp(k * deltaX);
                Aup = min(S * u, A);
				ku = round(log(Aup / runningMin) / deltaX, 0);
				Adown = min(S * d, A);
				kd = round(log(Adown / runningMin) / deltaX, 0);
                Vu = v(j+1+jshift, ku+kshift);
                Vd = v(j+jshift, kd+kshift);
                newV(j+jshift, k+kshift) = max(exp(-r * deltaT) * (p * Vu + (1-p) * Vd), K - A);
            end
        end
        v = newV;
    end
    
    optionValue = v(jshift, kshift);
    
end