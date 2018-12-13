clc; clear variables;

%% Midterm question 3
% Adam Frewin

% Constants
W = 6360; % weight
tsp = 1.2; % specific fuel consumption, lb/(lb*hr)
T_sl = 2040;
W_f = 2300; % fuel weight
S = 182; % wing area, ft
St = 42;
c = 5.47;
AR = 5.92; % aspect ratio


K = 0.08; % induced drag factor
e = 0.9; % oswald efficiency

Cd0 = 0.02; % drag
Cl_max = 1.4;
Cmac = 0.01;
Clow = 0.1;

a = 6.06; % lift curve slopes
at = 4.7; % tail lcs
ae = 0.4; % elevator evectiveness (slope of Cl vs delta_e)
l = 16; % distance between aerodynamic centers
h_cg = .10;
hac = .25;
it = -1*pi/180;
de_dalpha = 0.17;
e0 = 0;

rho = dens_imp(10000);
V_targ = 370; % ft/s

[alpha_req, delta_req] = static_stability_elevator(h_cg, hac, a, at, ae, S,...
    St, l, c, de_dalpha, Cmac, Clow, it, e0, rho, V_targ, W);