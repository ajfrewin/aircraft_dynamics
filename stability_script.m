clc; clear variables;

%% Constants
a = 6.5;
at = 5.3;
de_dalpha = .177;
e0 = 0;
it = -3 * pi / 180;
S = 290;
St = 70;
c = 6;
lt = 20;
Cmac = -.01;
Clow = 0;
h = 0.5;
h_ac = .25;
W = 35000;

rho = dens_imp(0);
% neutral point
h_np = neutral_point(a, at, S, St, lt, c, de_dalpha);

% trim AoA, no elevator
alpha_0 = static_stability(h, h_ac, a, at, S, St, lt, c, de_dalpha,...
    Cmac, Clow, it, e0, 0);

% Trim velocity
V_trim = sqrt(W /(.5*rho*S*alpha_0*a));

% elevator control
[alpha_trim, delta_trim] = static_stability_elevator(h, h_ac, a, at, 5,...
    S, St, lt, c, de_dalpha, Cmac, Clow, it, e0, rho, 500, W);