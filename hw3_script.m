clc; clear variables;

% Physical parameters
m = 13000/32.2; % mass
W = 13000; % weight
rho = dens_imp(40000); % air density at cruise altitude
S = 230; % wing area
c = 7; % chord length
thta0 = 0; % pitch at cruise
u0 = 677; % airspeed at cruise
Iy = 18800; % Moments of inertia
Ix = 28000;
Iz = 47000;
Izx = 1300;
b = 34; % wing span

% Longitudinal stability derivatives
CDu = .104;
CD1 = .0335;
CDalpha = .3;
CL1 = .41;
CLu = .4;
CLalpha = 5.84;
CLalphadot = 2.2;
CLq = 4.7;
Cm_u = .050;
Cm1 = 0;
Cmalpha = -.64;
Cmq = -15.5;
Cmalphadot = -6.7;

% Lateral stability derivatives
Clbeta = -.110;
Clp = -.45;
Clr = .16;
CYbeta = -.73;
CYp = 0;
CYr = .4;
Cnbeta = .127;
Cnp = -.008;
Cnr = -.2;

% Longitudinal Control Derivatives
CDdeltae = 0;
CTxu = -.07;
CTx1 = .0335;
CLdeltae = .46;
Cmdeltae = -1.24;
CmTu = -.003;
CmT1 = 0;

% Lateral Control Derivatives
Cldeltaa = .178;
Cldeltar = .019;
Cydeltaa = 0;
Cydeltar = .14;
CNdeltaa = -.02;
CNdeltar = -.074;

% State transition matricies
A_long = long_stab(CDu, CD1, CDalpha, CL1, CLu, CLalpha, CLalphadot,...
    CLq, Cm_u, Cm1,Cmalpha, Cmq, Cmalphadot, W, m, rho, S, c, thta0, u0, Iy);
fprintf('A_long\n');
disp(A_long);
A_lat = lat_stab(Clbeta, Clp, Clr, CYbeta, CYp, CYr, Cnbeta, Cnp, Cnr,...
    m, rho, S, b, thta0, u0, Ix, Iz, Izx);
fprintf('A_lat\n');
disp(A_lat);

% Control Matricies
Czalpha = -(CLalpha + CD1);
Zwdot = m-1/2*rho*u0*S*Czalpha/A_long(2,2);
Mwdot = 1/4*rho*c^2*S*Cmalphadot;
B_long = long_cntrl(CDdeltae, CTxu, CTx1, CLdeltae, Cmdeltae, CmTu,...
    CmT1, rho, u0, S, m, c, Zwdot, Mwdot, Iy);
fprintf('B_long\n');
disp(B_long);
B_lat = lat_cntrl(Cydeltaa, Cydeltar, Cldeltaa, Cldeltar, CNdeltaa, CNdeltar,...
    rho, u0, S, b, m, Ix, Iz, Izx);
fprintf('B_lat\n');
disp(B_lat);


long_modes = eig(A_long);
lat_modes = eig(A_lat);
fprintf('Long modes\n');
disp(long_modes);
fprintf('Lat modes\n');
disp(lat_modes);

%% Longitudinal "stick fixed" response
[tsim_long, x_long] = ode45(@(t,x) A_long*x, [0 400], [10 0 0 0]');
figure('Name','Longitudinal "Stick Fixed" response');
subplot(4,1,1);
plot(tsim_long, x_long(:,1)); grid(); ylabel('\Delta u ft/s');

subplot(4,1,2);
plot(tsim_long, x_long(:,2)); grid(); ylabel('\Delta w ft/s');

subplot(4,1,3);
plot(tsim_long, x_long(:,3)*180/pi); grid(); ylabel('\Delta q ^\circ /s');

subplot(4,1,4);
plot(tsim_long, x_long(:,4)*180/pi); grid(); ylabel('\Delta \theta ^\circ');

%% Lateral "stick fixed" response
figure('Name','Lateral "stick fixed" response')
[tsim_lat, xsim_lat] = ode45(@(t,x) A_lat*x, [0 400], [0 1*pi/180 0 0]');
subplot(4,1,1);
plot(tsim_lat, xsim_lat(:,1)); grid(); ylabel('v ft/s');

subplot(4,1,2);
plot(tsim_lat, xsim_lat(:,2)); grid(); ylabel('p ^\circ /s');

subplot(4,1,3);
plot(tsim_lat, xsim_lat(:,3)); grid(); ylabel('r ^\circ /s');

subplot(4,1,4);
plot(tsim_lat, xsim_lat(:,4)); grid(); ylabel('phi ^\circ');

%% Forced response
long_sys_elev = ss(A_long, B_long(:,1), eye(4), zeros(4,1));
long_sys_thrust = ss(A_long, B_long(:,2), eye(4), zeros(4,1));

figure('Name', "Longitudinal elevator impulse response");
impulse(long_sys_elev); grid();
figure('Name', "Longitudinal elevator step response");
step(long_sys_elev); grid();
figure('Name', "Longitudinal thrust impulse response");
impulse(long_sys_thrust); grid();
figure('Name', "Longitudinal thrust step response");
step(long_sys_thrust); grid();

lat_sys = ss(A_lat, B_lat, eye(4), zeros(4,2));
figure('Name', "Lateral impulse response");
impulse(lat_sys);
figure('Name', "Lateral step response");
step(lat_sys);

% hold at 5 degrees
Q = .01 * eye(4);
R = 10*eye(2);
K_lat_lqr = lqr(A_lat, B_lat, Q, R, zeros(4,2));
disp(eig(A_lat - B_lat*K_lat_lqr));
[tsim_lat_cntrl, xsim_lat_cntrl] = ode45(@(t,x) (A_lat - B_lat*K_lat_lqr)*...
    (x - [0; 0; 0; 5*pi/180]), [0 400], [0 1*pi/180 0 0]');

figure('Name','5 degree roll angle control')
subplot(4,1,1);
plot(tsim_lat_cntrl, xsim_lat_cntrl(:,1)); grid(); ylabel('v ft/s');

subplot(4,1,2);
plot(tsim_lat_cntrl, xsim_lat_cntrl(:,2)); grid(); ylabel('p ^\circ /s');

subplot(4,1,3);
plot(tsim_lat_cntrl, xsim_lat_cntrl(:,3)); grid(); ylabel('r ^\circ /s');

subplot(4,1,4);
plot(tsim_lat_cntrl, xsim_lat_cntrl(:,4)*180/pi); grid(); ylabel('phi ^\circ');

%% Altitude controller
thta0 = 0;
u0 = 677;
A_long_alt = [A_long zeros(4,1); -sin(thta0) cos(thta0) 0 -u0*cos(thta0) 0];
    
K_alt_hold = lqr(A_long_alt, [B_long; 0 0], .01*eye(5), 10*eye(2), zeros(5,2));

