% Author: Adam Frewin Nov 2018
% Description: script to determine absolute and service ceilings of an
% aircraft
% Data for this script is from the Gulfstream IV
clc; clear variables;

%% constants
W = 73000; % weight
S = 950; % wing area
Cd0 = 0.015; % Constand drag coefficient
K = 0.08; % induced drag factor
AR = 5.92; % aspect ratio
T_sl = 2*13850; % Thrust at sea-level
rho0 = 2.3760e-3; % sea level density

L_D_max = sqrt(1/(4*Cd0*K)); % Max L / D
hs = 10000:100:100000; % Creating for graph visualizations, span of altitudes
rc_max = zeros(1,length(hs));

% Calculating RC_max as function of altitude
for i=1:length(hs)
    rhos = dens_imp(hs(i));
    t = T_sl * rhos/rho0;
    z = 1 + sqrt(1 + 3 / (L_D_max^2 * (t/W)^2));
    rc_max(i) = sqrt((W/S * z)/(3*rhos*Cd0)) * (t/W)^(3/2) * ...
    (1 - z/6 - 3/(2*(t/W)^2 * L_D_max^2*z));
end
combined = [hs; rc_max]; % for ease of matching altitude to RoC



%% Iterative method for finding absolute cieling

eps = 5; % convergence threshold
dh = 5; % change in height, need to compute 2 points per iteration
h = 10000; % initial height
RC_max1 = eps+1; % initializing RC_max, value will be reassigned

while abs(RC_max1) > eps
    
    % Newton-Raphson method requires function derivative for iteration,
    % in this case derivative is mathematically intense. RC_max vs h is 
    % quite linear in this regime, as seen in the graph, so we can 
    % use a numeric approximation for the derivative at the given point
    
    rho1 = dens_imp(h);
    rho2 = dens_imp(h+dh);
    T1 = T_sl * rho1/rho0;
    T2 = T_sl * rho2/rho0;
    
    z1 = 1 + sqrt(1 + 3 / (L_D_max^2 * (T1/W)^2));
    RC_max1 = sqrt((W/S * z1)/(3*rho1*Cd0)) * (T1/W)^(3/2) * ...
    (1 - z1/6 - 3/(2*(T1/W)^2 * L_D_max^2*z1));

    z2 = 1 + sqrt(1 + 3 / (L_D_max^2 * (T2/W)^2));
    RC_max2 = sqrt((W/S * z2)/(3*rho2*Cd0)) * (T2/W)^(3/2) * ...
    (1 - z2/6 - 3/(2*(T2/W)^2 * L_D_max^2*z2));
    
    % numeric estimation of derivative
    RC_prime = (RC_max2 - RC_max1) / dh;
    
    delta = -RC_max1/RC_prime;
    h = h + delta;
    disp(h);
end

% visual of rc_max as function of altitude
plot(hs, rc_max, 'Linewidth', 2)
hold on
plot(hs, zeros(1,length(hs)),'r--');
plot([h h], [min(rc_max)-2 max(rc_max)+2], 'r--')
ylim([min(rc_max) max(rc_max)]);
grid();
xlabel('altitude (ft)');
ylabel('RC_{max} (ft)');

%% Re-initialize variables for service ceiling
h = 10000; % initial height
RC_max1 = eps+1; % initializing RC_max, value will be reassigned
while abs(RC_max1) > eps
    
    % Same method as above,
    % just solves for when RC_max = 1.667 ft/s = 100 ft/min
    
    rho1 = dens_imp(h);
    rho2 = dens_imp(h+dh);
    T1 = T_sl * rho1/rho0;
    T2 = T_sl * rho2/rho0;
    
    z1 = 1 + sqrt(1 + 3 / (L_D_max^2 * (T1/W)^2));
    RC_max1 = sqrt((W/S * z1)/(3*rho1*Cd0)) * (T1/W)^(3/2) * ...
    (1 - z1/6 - 3/(2*(T1/W)^2 * L_D_max^2*z1))-1.667;

    z2 = 1 + sqrt(1 + 3 / (L_D_max^2 * (T2/W)^2));
    RC_max2 = sqrt((W/S * z2)/(3*rho2*Cd0)) * (T2/W)^(3/2) * ...
    (1 - z2/6 - 3/(2*(T2/W)^2 * L_D_max^2*z2))-1.667;
    
    % numeric estimation of derivative
    RC_prime = (RC_max2 - RC_max1) / dh;
    
    delta = -RC_max1/RC_prime;
    h = h + delta;
    disp(h);
end


