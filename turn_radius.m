clc; clear variables;

W = 73000; % weight
S = 950; % wing area
Cd0 = 0.015; 
K = 0.08; % induced drag factor
h = 0; % altitude
e = 0.9; % oswald efficiency
AR = 5.92; % aspect ratio
T_sl = 1; % Thrust at sea-level
g = 32.17;
rho = dens_imp(h);

T_W = .178;
L_D_max = 14.46;
W_S = 76.84;
CL_max  = 1.2;

V = 200:100:800;
n_max_alph = .5*rho*V.^2 * CL_max / (W/S);
n_max_T = sqrt(.5*rho*V.^2/(K*W/S).* (T_W - .5*rho*V.^2*Cd0/(W/S)));
n_max_S = 2.5*ones(1,length(V));
n_lim = min(n_max_alph, n_max_T);
n_lim = min(n_lim, n_max_S);

R_min = V.^2.*(g*sqrt(n_lim.^2 - 1)).^-1;


fprintf('   V            N alpha\n');
disp('-----------------------');
count = 200;
diff = V(2) - V(1);
for i=1:length(V)
    disp(count + " ft/s : " + n_max_alph(i));
    count = count + diff;
end

fprintf('   V           N T\n');
disp('-----------------------');
count = 200;
diff = V(2) - V(1);
for i=1:length(V)
    disp(count + " ft/s : " + n_max_T(i) + ' ft/s');
    count = count + diff;
end

plot(V, n_max_alph);
hold on
plot(V, n_max_T);