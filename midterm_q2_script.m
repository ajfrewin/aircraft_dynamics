clc; clear variables;

%% Midterm Question 2
%  Adam Frewin

% Constants
W = 6360; % weight
ct = 1.2; % specific fuel consumption, lb/(lb*hr)
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

a = 6.06;
at = 4.7;
l = 16;
h_cg = .10;
it = -1*pi/180;
de_dalpha = 0.17;

hs = [5000 10000]; % altitude (ft)
% atmosphere conditions
fig_clcd = figure();
for i=1:length(hs)
    h = hs(i);
    tau = temp_ft(h);
    rho = dens_imp(h);

    % velocity vector
    V = 100:10:600;

    % coeff of lift
    Cl = W * (.5*rho*V.^2*S).^-1;

    % drag coeff
    Cd = Cd0 + K*Cl.^2;

    % glide coeff
    Cl_Cd = Cl.*Cd.^-1;


    % Max Values
    L_D_max = sqrt(1/(4*Cd0*K));
    V_LD_max = sqrt(W / (.5*rho*S*sqrt(Cd0/K)));

    % Display important values
    disp('Altitude = ' + string(h));
    disp('W/S = ' + string(W/S));
    disp('L/D max = ' + string(L_D_max));
    disp('V(L/D max) = ' + string(V_LD_max));
    fprintf('\n');
    
    % Plots
    figure(fig_clcd);
    plot(V,Cl_Cd);
    hold on
end
%% Now make them plots pretty
figure(fig_clcd);
title('L / D vs V');
xlabel('V (ft/s)');
ylabel('L / D');
legend('h = 5,000 ft', 'h = 10,000 ft', 'Location', 'southeast');
grid();

%% RC max calc
rho0 = dens_imp(0);
L_D_max = sqrt(1/(4*Cd0*K)); % Max L / D
hs = 10000:100:100000; % Creating for graph visualizations
rc_max = zeros(1,length(hs));
% Calculating RC_max as function of altitude
for i=1:length(hs)
    rhos = dens_imp(hs(i));
    t = T_sl * rhos/rho0;
    z = 1 + sqrt(1 + 3 / (L_D_max^2 * (t/W)^2));
    rc_max(i) = sqrt((W/S * z)/(3*rhos*Cd0)) * (t/W)^(3/2) * ...
    (1 - z/6 - 3/(2*(t/W)^2 * L_D_max^2*z));
end
combined = [hs; rc_max]; % for ease of matching height to value

%% Iterative method for finding absolute cieling

eps = 1; % convergence threshold
dh = 5; % change in height, need to compute 2 points per iteration
ha = 10000; % initial height
RC_max1 = eps+1; % initializing RC_max, value will be reassigned
fprintf('Calculating absolute ceiling:\n');
while abs(RC_max1) > eps
    
    % Newton-Raphson method requires function derivative for iteration,
    % in this case derivative is mathematically intense. RC max is 
    % quite linear in this regime, as seen in the graph, so we can 
    % use a numeric approximation for the derivative at the given point
    
    rho1 = dens_imp(ha);
    rho2 = dens_imp(ha+dh);
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
    ha = ha + delta;
    disp(ha);
end

% visual of rc_max as function of altitude
figure();
plot(hs, rc_max, 'Linewidth', 2)
hold on
plot(hs, zeros(1,length(hs)),'r--');
plot([ha ha], [min(rc_max)-2 max(rc_max)+2], 'r--')
ylim([min(rc_max) max(rc_max)]);
grid();
xlabel('altitude (ft)');
ylabel('RC_{max} (ft)');

%% Re-initialize variables for service ceiling
hserv = 30000; % initial height
RC_max1 = eps+1; % initializing RC_max, value will be reassigned
fprintf('Calculating service ceiling: \n');
while abs(RC_max1) > eps
    
    % Similar to above section, just solves for when RC_max = 1.667 ft/s =
    % 100 ft/min
    
    rho1 = dens_imp(hserv);
    rho2 = dens_imp(hserv+dh);
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
    hserv = hserv + delta;
    disp(hserv);
end

%% Minimum turn radius calc
T_W = T_sl/W; % thrust to weight
L_D_max = sqrt(1/(4*Cd0*K)); % max L/D
g = 32.17;
V_targ = 200;
n_max_alph = .5*rho*V_targ^2 * Cl_max / (W/S); % lifting limit on load factor
n_max_T = sqrt(.5*rho*V_targ^2/(K*W/S).* (T_W - .5*rho*V_targ^2*Cd0/(W/S)));
% limit on load factor from available thrust
% no structural limit
n_lim = min(n_max_alph, n_max_T);

R_min = V_targ^2/(g*sqrt(n_lim^2 - 1));
disp('R_min = ' + string(R_min));

%% Range calculation
% velocity vector
rho_2 = dens_imp(10000);
V = 100:10:600;
% coeff of lift
Cl = W * (.5*rho_2*V.^2*S).^-1;
% drag coeff
Cd = Cd0 + K*Cl.^2;
% glide coeff
Cl_Cd = Cl.*Cd.^-1;

V_times_glide = V.*Cl_Cd;
[val, idx] = max(V_times_glide);
maximizing_V = V(idx);
maximizing_ClCd = Cl_Cd(idx);
disp('Range-maximizing Velocity = ' + string(maximizing_V) + ' ft/s');

Rng = val / ct * log(W/(W-W_f)) * 3600;
disp('Range = ' + string(Rng/5280) + ' miles')



