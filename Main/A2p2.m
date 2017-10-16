
explicitIIIFDM(9.8, 9, 0.001, 0.25, 0.01, 0.15, 25, 0.05)

S0 = 9.8;
X = 9;
T = 0.25;
r = 0.001;
q = 0.01;
sigma = 0.15;
N = 25;
ds = 0.05;

Smax = 4 * X;
I = floor(Smax / ds);

N = T * (sigma ^ 2 * I ^ 2 + r * I);

explicitIIIFDM(9.8, 9, 0.001, 0.25, 0.01, 0.15, 2916, 0.05)


explicitIIIFDMAmericanCall(9.8, 9, 0.001, 0.25, 0.01, 0.15, 2900, 0.05)