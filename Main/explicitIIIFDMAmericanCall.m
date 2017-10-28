%Group members: Chen Penghao, Wang Zexin
%Group number: G01

function [optionValue] = explicitIIIFDMAmericanCall(S0, X, r, T, q, sigma, N, ds)
    % Initialize parameters
    Smax = 4 * X;
    dt = T / N;
    I = round(Smax / ds);
    
    % Initialize discretized lattice
    VGrid = zeros(I+1, N+1);
    VGrid(1, :) = 0;
    VGrid(I+1, :) = (Smax-X) * exp(-r*(T:-dt:0));
    
    % Initialize terminal payoff
    VGrid(:, N+1) = max((0:I) * ds - X, 0);
    
    % Initialize explicit scheme III coefficients
    i = (1 : (I-1))';
    isq = i.^2;
    c = (0.5 * sigma ^ 2 * isq + (r-q) * i) * dt / (1 + r * dt);
    b = (1 - sigma ^ 2 * isq * dt - (r-q) * i * dt) / (1 + r * dt);
    a = (0.5 * sigma ^ 2 * isq) * dt / (1 + r * dt);
    
    % Verify positivity condition
    len01 = length(a);
    len02 = length(find(a < 0));
    disp(['Coeff a, Of ', num2str(len01), ' elements, ', num2str(len02), ' violated the positivity condition.']);
    len01 = length(b);
    len02 = length(find(b < 0));
    disp(['Coeff b, Of ', num2str(len01), ' elements, ', num2str(len02), ' violated the positivity condition.']);
    len01 = length(c);
    len02 = length(find(c < 0));
    disp(['Coeff c, Of ', num2str(len01), ' elements, ', num2str(len02), ' violated the positivity condition.']);
    
    % Backward iteration
    ishift = 1;
    for n = N : -1 : 1
        row = a .* VGrid(i - 1 + ishift, n + 1) + b .* VGrid(i + ishift, n + 1) + c .* VGrid(i + 1 + ishift, n + 1);
        payoff = (i * ds - X) .* exp(-r * n * dt);
        row = max(row, payoff);
        VGrid(i + ishift, n) = row;
    end
    
    % Obtain the option value
    optionValue = VGrid(round(S0 / ds) + ishift, 1);

end

function optionPrice =EuropeanVanillaCall(S, K, r, tau, sigma, q)
	d1 = (log(S/K) + (r-q+sigma*sigma/2) * tau) / sigma / sqrt(tau);
	d2 = d1 - sigma * sqrt(tau);
	optionPrice = S*normcdf(d1) - normcdf(d2)*X*exp(-r*tau);
end