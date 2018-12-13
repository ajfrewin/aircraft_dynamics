clc; clear variables;

%% Constants
W = 13000; % weight
alt = 30000; % altitude
h = 0.25; % cg loaction
Cmac = -.01;
c = 5.86; % mean chord
Cd0_t = .02; % tail drag
Kt = .062; % tail induced drag coeff
it = -3*pi/180; % tail incidence
ae = 0.5; % elevator effectiveness
span = 53 + 4/12; % wing span
span_t = 20.83; % tail span
S = 294; % wing area
St = 70.68; % tail area
sweep = 4*pi/180; % wing sweep
sweep_t = 23*pi/180; % tail sweep
ltv = 8.8;
lth = 45.1;
lt = sqrt(ltv^2+lth^2);
AR = span^2/S; % wing AR
AR_t = span_t^2/St; % tail AR
Clow = 0; 
e0 = 0;
Cmac_t = 0;
Vh = St/S *lt/c;
taper = 2.5/8.2;
rho = dens_imp(alt);
hac = .25;
SoS = 678.1*5280/3600; % speed of sound at 30,000 ft
h = .25;
Mach = 0.05:0.01:.98;
V = 675;
%{
a_w = lift_curve(AR, Mach, sweep);
a_t = lift_curve(AR_t, Mach, sweep_t);
de_dalpha = 4.44*sqrt(1-Mach.^2).*((1/AR - 1/(1 + AR^1.7))*(10-3*taper)/7.*...
    (1 - ltv/span)/((2*lt/span)^.33)*sqrt(cos(sweep))).^1.19;
%}
a_w = lift_curve(AR, 0, sweep);
a_t = lift_curve(AR_t, 0, sweep_t);
de_dalpha = 4.44*sqrt(1).*((1/AR - 1/(1 + AR^1.7))*(10-3*taper)/7.*...
    (1 - ltv/span)/((2*lt/span)^.33)*sqrt(cos(sweep))).^1.19;
a_bar = a_w + a_t.*(1 - de_dalpha).*St/S;
h_np = hac + a_t./a_w.*(1 - de_dalpha).*Vh;


theta = atan(ltv/lth);
alphas = 0:.01:5;
alphas = alphas*pi/180;
alpha_t = alphas - de_dalpha.*alphas + it;
deltas = [-15 -10 -5 0 5 10 15];
figure();
for i=1:length(deltas)
    delta = deltas(i)*pi/180;
    Cm = Cmac - St/S*lt*a_t.*alpha_t - St/S*lt.*(Cd0_t + Kt.*...
        (a_t.*alpha_t).^2).*(alpha_t - theta) - ae*Vh*delta;
    [val, idx] = min(abs(Cm));
    fprintf('Delta_e = ' + string(delta*180/pi) + ' degrees, '+...
        'Trim alpha = ' + string(alphas(idx)*180/pi) + ' degrees\n');
    plot(alphas*180/pi, Cm);
    hold on
end
fprintf('\n');
grid();
plot([0 5], [0 0], 'k--');
plot([0 0], [-2 3], 'k--');
legend(string(deltas) + ' degrees');
xlabel('\alpha');
ylabel('C_m');

a_w = lift_curve(AR, V/SoS, sweep);
a_t = lift_curve(AR_t, V/SoS, sweep_t);
de_dalpha = 4.44*sqrt(1-(V/SoS)^2).*((1/AR - 1/(1 + AR^1.7))*(10-3*taper)/7.*...
    (1 - ltv/span)/((2*lt/span)^.33)*sqrt(cos(sweep))).^1.19;
[alpha, delta] = static_stability_elevator(h, hac, a_w, a_t, ae,...
    S, St, lt, c, de_dalpha, Cmac, Clow, it, e0, rho, V, W);