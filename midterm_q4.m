clc; clear variables;

%% Midterm question 4
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

A = 2;
V0 = 160;
y0 = 30;
z0 = -55;
b = [0 1.38578e-3 -.114569 3.56385 -24.8365 42.8718 30];
c = [0 -4.03145e-4 3.29785e-2 -.976344 4.54341 8.68104 -55];

ts = 0:.01:4.64466;
x_t = V0*ts - A*ts.^2;
y_t = b(2) * ts.^5 + b(3)*ts.^4 + b(4)*ts.^3 + b(5)*ts.^2 + b(6)*ts + b(7);
z_t = c(2) * ts.^5 + c(3)*ts.^4 + c(4)*ts.^3 + c(5)*ts.^2 + c(6)*ts + c(7);

plot3(x_t, y_t, z_t);
hold on
scatter3(x_t(1), y_t(1), z_t(1), 20, 'r');
grid();
xlabel('x');
ylabel('y');
zlabel('z');