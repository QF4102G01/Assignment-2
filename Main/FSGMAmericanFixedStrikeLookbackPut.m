function [optionValue] = FSGMAmericanFixedStrikeLookbackPut(runningTime, tau, S0, sigma, q, runningMin, r, K, N)

    deltaT = tau / N;
    deltaX = sigma * sqrt(deltaT);
    u = exp(deltaX);
    d = exp(-deltaX);
    p = (exp(r * deltaT) - d) / (u-d);
    jshift = 1;
    kshift = N + 1;
    runningMin = min(runningMin, S0);

    v = zeros(N+1, 2 * N + 1);
    for j = 0 : N
        for k = -N : N
            A = S0 * exp(k * deltaX);
            v(j+jshift, k+kshift) = max(K - A, 0);
        end
    end
    
    for n = N-1 : -1 : 0
            newV = zeros(N+1, 2 * N + 1);
            aS = zeros(N+1, 2 * N + 1);
            for j = 0 : n
                for k = -n : n
                    S = S0 * exp((2 * j - n) * deltaX);
                    A = S0 * exp(k * deltaX);

                    aS(j+jshift, k+kshift) = S;

                    ku = k;
                    kd = min(k, 2 * j - n);
                    Vu = v(j+1+jshift, ku+kshift);
                    Vd = v(j+jshift, kd+kshift);
                    newV(j+jshift, k+kshift) = max(exp(-r * deltaT) * (p * Vu + (1-p) * Vd), K - A);
                end
            end
        v = newV;
        aS;
    end
    
    k = floor(log(runningMin / S0) / deltaX);
    Au = S0 * exp((k+1) * deltaX);
    Ad = S0 * exp(k * deltaX);
    alpha = (Au - runningMin) / (Au - Ad);
    
    optionValue = max(alpha * v(jshift, kshift + k) + (1 - alpha) * v(jshift, kshift + k + 1), K - runningMin);
    
end