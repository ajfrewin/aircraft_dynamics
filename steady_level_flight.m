clc; clear variables;

%% Determines Glide coefficient and thrust required
%  for steady level flight
%  Adam Frewin and Mason Handy October 2018

% Constants
W = 73000; % weight
S = 950; % wing area
Cd0 = 0.015; 
K = 0.08; % induced drag factor
e = 0.9; % oswald efficiency
AR = 5.92; % aspect ratio

hs = [10000 30000 50000]; % altitude (ft)
% atmosphere conditions
fig_clcd = figure();
fig_tr = figure();
for i=1:3
    h = hs(i);
    tau = temp_ft(h);
    rho = dens_imp(h);

    % velocity vector
    V = 200:10:1000;

    % coeff of lift
    Cl = W * (.5*rho*V.^2*S).^-1;

    % drag coeff
    Cd = Cd0 + K*Cl.^2;

    % glide coeff
    Cl_Cd = Cl.*Cd.^-1;

    %Thrust required
    Tr = W * Cl_Cd.^-1;

    % Max Values
    L_D_max = sqrt(1/(4*Cd0*K));
    V_LD_max = sqrt(W / (.5*rho*S*sqrt(Cd0/K)));

    % Display important values
    disp('Altitude = ' + string(h));
    disp('W/S = ' + string(W/S));
    disp('L/D max = ' + string(L_D_max));
    disp('V(L/D max) = ' + string(V_LD_max));
    fprintf('\n');
    
    Tr_min = W / L_D_max;
    
    % Plots
    figure(fig_clcd);
    plot(V,Cl_Cd);
    hold on
    figure(fig_tr);
    plot(V, Tr);
    hold on
end
%% Now make them plots pretty
figure(fig_clcd);
title('L / D vs V');
xlabel('V (ft/s)');
ylabel('L / D');
legend('h = 10,000 ft', 'h = 30,000 ft', 'h = 50,000 ft',...
    'Location', 'southeast');
grid();

figure(fig_tr);
title('Required Thrust vs V');
xlabel('V (ft/s)');
ylabel('T_{req} (lbf)');
legend('h = 10,000 ft', 'h = 30,000 ft', 'h = 50,000 ft',...
    'Location', 'northeast');
grid();
%% Text outputs
disp('Cl / Cd');
count = 200;
diff = V(2) - V(1);
for i=1:length(Cl)
    %disp(count + " ft/s : " + Cl_Cd(i));
    count = count + diff;
end
disp('Thrust Req');
count = 200;
for i=1:length(Cl)
    %disp(count + " ft/s : " + Tr(i) + ' lbf');
    count = count + diff;
end










