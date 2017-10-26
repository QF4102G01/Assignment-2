%Group members: Chen Penghao, Wang Zexin
%Group number: G01

function [optionValue] = explicitIIIFDMAmericanCall(S0, X, r, T, q, sigma, N, ds)

    Smax = 4 * X;
    dt = T / N;
    I = round(Smax / ds);
    VGrid = zeros(I+1, N+1);
    
    VGrid(1, :) = X * exp(-r * (T - 0 : dt : T));
    VGrid(I+1, :) = 0;
    
    VGrid(:, N+1) = max(X - (0:I) * ds, 0);
    
    i = (1 : (I-1))';
    isq = i.^2;
    c = (0.5 * sigma ^ 2 * isq + (r-q) * i) * dt / (1 + r * dt);
    b = (1 - sigma ^ 2 * isq * dt - (r-q) * i * dt) / (1 + r * dt);
    a = (0.5 * sigma ^ 2 * isq) * dt / (1 + r * dt);
    
    len01 = length(a);
    len02 = length(find(a < 0));
    disp(['Coeff a, Of ', num2str(len01), ' elements, ', num2str(len02), ' violated the positivity condition.']);
    len01 = length(b);
    len02 = length(find(b < 0));
    disp(['Coeff b, Of ', num2str(len01), ' elements, ', num2str(len02), ' violated the positivity condition.']);
    len01 = length(c);
    len02 = length(find(c < 0));
    disp(['Coeff c, Of ', num2str(len01), ' elements, ', num2str(len02), ' violated the positivity condition.']);
    
    ishift = 1;
    for n = N : -1 : 1
        row = a .* VGrid(i - 1 + ishift, n + 1) + b .* VGrid(i + ishift, n + 1) + c .* VGrid(i + 1 + ishift, n + 1);
        payoff = max(i * ds - X, 0);
        row = max(row, payoff);
        VGrid(i + ishift, n) = row;
    end
    
    optionValue = VGrid(round(S0 / ds) + ishift, 1);

end