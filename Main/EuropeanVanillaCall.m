%Group members: Chen Penghao, Wang Zexin
%Group number: G01

function [c] = EuropeanVanillaCall(S0, q, X, tau, r, sigma)
  d1 = (log(S0./X) + (r - q + sigma ^ 2 / 2) * tau) / (sigma * tau ^ 0.5);
  d2 = d1 - sigma * tau ^ 0.5;
  c = S0.* exp(-q * tau).* normcdf(d1) - X * exp(-r * tau) * normcdf(d2);
end
