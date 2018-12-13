clear; close all; clc;

temp_sl = 288.16;       % Temperature at sea level  K
press_sl = 1.01325;     % Pressure at sea level     N/m^2
density_sl = 1.225;     % Density at sea level      kg/m^3
alt = 9144;             % altitude                  ft
a1 = -6.5*10^-3;        % Temperature gradient      K/m
g = 9.8;                % gravity                   g
R = 287;                % Specific gas constant     J/(Kg*K)
temp_h = temp_sl + (a1*alt);                             %Temperature at h K
press_h = press_sl*(temp_h/temp_sl).^(-g/(a1*R));        %Pressure at h N/m^2
density_h = density_sl*(temp_h/temp_sl).^(-1-g/(a1*R));  %Density at h kg/m^3

density_h_imp = density_h * 0.0019403203;

% Values   
W = 13000;                      % Weight                    lb
cruise = 30000;                 % Cruising altitude         ft
CG = 0.25;                      % CG Location               ft
C_LOW = -0.01;                  % Zero lift pitching moment of the wing
c = 5.86;                       % Mean wing chord length    ft
C_D0 = 0.02;                    % Horizontal tail drag coefficient
K = 0.062;                      % Horiz tail drag coeff
i = -3;                         % Horiz. tail incidence angle% degrees
ee = 0.5;                       % elevator effectiveness    

% Derived values 
v = 675;                        % Airspeed                  ft/s
sound = sqrt(1.4*1716*411.7032);
AR = 26.2/5.86;                 % Aspect Ratio              ft   2?
AR_t = 20.83/((4.7+2)/2);         % Aspect Ratio              ft
mach = v/sound;                 % Mach number 
% geometric values
b = 53.33;                      % Wing Span                 f
St = .5*(2+4.7)*.5*20.83;       % Surface Area of Tail      ft^2
S = 26.2*5.86;                  % Wing area                 ft^2
lbd = 2.5/8.2;                  % Taper Ratio of Wing
lbdt = 2/4.7;                   % Taper Ratio of Tail 
Lambda = deg2rad(4);            % Sweep Angle of wing       rad
Lambda_t = deg2rad(23);            % Sweep Angle of tail
MACt = 4.7*2/3*((1+lbdt+lbdt^2) / (1+lbdt));   % Mean Aerodynamic Chord of Tail 
h_ac = 0.25*5.86;               % Aerodynamic Center Distance for Wing    ft
h_ac_t = 0.25*MACt;             % Aerodynamic Center Distance for Tail    ft
l_th = 23.3-h_ac+h_ac_t;        % Horizontal distance between wing AC and tail AC 
l_tv = 11.7;                    % vertical distance between wing AC and tail AC 
l_t = sqrt(l_th^2 + l_tv^2);    % Distance from wing AC to tail AC 
Vh_h = l_th*St/c/S;             % Horizontal Tail volume ratio
Vh_t = l_tv*St/c/S;             % Vertical Tail volume ratio
Vh = l_t*St/c/S;                % Real Tail colume ratio 
%% Part A

k = 1+(8.2-2.3*lbd-AR*(0.22-.153*lbd))/100;        % constant for lift curve slope 
k_t = 1+((AR_t)*(1.87-2.33*10^(-4)*Lambda_t))/100;

lift_curve_slope = (2*pi*AR) / (2+ sqrt(((AR^2*(1-mach^2))/k^2)*(1+(tan(Lambda)^2)/(1-mach^2))+4));
lift_curve_slope_t = (2*pi*AR_t) / (2+ sqrt(((AR_t^2*(1-mach^2))/k_t^2)*(1+(tan(Lambda_t)^2)/(1-mach^2))+4));

d_eps_d_alph = 4.44*sqrt(1-mach^2)*((1/AR - 1/(1+AR^(1.7)))*((10-3*lbd)/7)*((1-l_tv/b)/(2*l_th/b)^0.33)*sqrt(cos(Lambda)))^1.19;
a_bar = lift_curve_slope + lift_curve_slope_t*(1-d_eps_d_alph)*St/S;

h_np = 0.25 + (lift_curve_slope_t/a_bar)*(1-d_eps_d_alph)*Vh;

SS = h_np - 0.25;
SS2 = SS*c;

syms C_mt(alpha)
C_mt(alpha) = Vh_h*((lift_curve_slope_t*(alpha-(d_eps_d_alph*alpha)+i)) + (C_D0+K*(lift_curve_slope_t*(alpha-(d_eps_d_alph*alpha)+i))^2)*(alpha-(d_eps_d_alph*alpha))) - Vh_t*((C_D0+K*(lift_curve_slope_t*(alpha-(d_eps_d_alph*alpha)+i))^2 - lift_curve_slope_t*(alpha-(d_eps_d_alph*alpha)+i))*(alpha-(d_eps_d_alph*alpha)));

dC_mt_dalph = diff(C_mt,alpha);


