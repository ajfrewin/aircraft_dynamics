% Author: Adam Frewin, Oct 2018
% Description: Script that calculates various values concerning steady
% climb and plots rate of climb vs velocity, angle of climb vs velocity
% Data taken from Gulfstream IV
clc; clear variables;

%% constants
W = 73000; % weight
S = 950; % wing area
Cd0 = 0.015; 
K = 0.08; % induced drag factor
h = 30000; % altitude
e = 0.9; % oswald efficiency
AR = 5.92; % aspect ratio
T_sl = 0.1; % Thrust at sea-level
    
% altitude dependent values
rho = dens_imp(h);
T = rho / 2.3769e-3 * T_sl;

% Max parameters
L_D_max = sqrt(1/(4*Cd0*K));
V_LD_max = sqrt(W / (.5*rho*S*sqrt(Cd0/K)));

% Max Rate of Climb calculation
z = 1 + sqrt(1 + 3 / (L_D_max^2 * (T/W)^2));
RC_max = sqrt((W/S * z)/(3*rho*Cd0)) * (T/W)^(3/2) * ...
    (1 - z/6 - 3/(2*(T/W)^2 * L_D_max^2*z));

% Note: due to the discretization, T cannot equal zero
disp('T/W = ' + string(T/W));
disp('W/S = ' + string(W/S));
disp('L/D max = ' + string(L_D_max));
disp('V(L/D max) = ' + string(V_LD_max));
disp('RC_max = ' + string(RC_max));
fprintf('\n');

% Vector span of velocities
V = 200:10:1000;

% coeff of lift
Cl = W * (.5*rho*V.^2*S).^-1;

% drag coeff
Cd = Cd0 + K*Cl.^2;

% glide coeff
Cl_Cd = Cl.*Cd.^-1;
    
% Small Angle approximation for RC
RC = V.*(T / W - 1/2*rho*V.^2*Cd0 / (W/S) - ...
    W / S * K *(.5*rho*V.^2).^-1);
gma = asin(RC./V); % Angle of climb

%% Value Printing
fprintf('   V           RC\n');
disp('-----------------------');
count = 200;
diff = V(2) - V(1);
for i=1:length(V)
    disp(count + " ft/s : " + RC(i) + ' ft/s');
    count = count + diff;
end

%% Plots
% RoC vs V
plot(V,RC);
hold on
plot(V, zeros(length(V)) + RC_max, 'r--');
plot([480 480], [-100 RC_max + 20], 'r--');
grid();
xlabel('V (ft/s)');
ylabel('Rate of climb (ft/s)');

figure();
plot(V, gma*180/pi);
grid();
ylabel('\gamma (degrees)');
xlabel('V (ft/s)');
