function A_long = long_stab(CDu, CD1, CDalpha, CL1, CLu, CLalpha, CLalphadot,...
    CLq, Cm_u, Cm1, Cmalpha, Cmq, Cmalphadot, W, m, rho, S, c, thta0, u0, Iy)
% Author: Adam Frewin 2018
% Description: Computes longitudinal stability matrix of a standard fixed-wing
% aircraft, according to the longitudinal linear model

g = 32.26;  % gravity acceleration, ft/s^2

% Non-dimensional stability derivatives
Cxu = -(CDu + 2*CD1);
Cxalpha = -(CDalpha - CL1);
Cxalphadot = 0;
Cxq = 0;
Czu = -(CLu + 2*CL1);
Czalpha=-(CLalpha + CD1);
Czalphadot = -CLalphadot;
Czq = -CLq;
Cmu = Cm_u + 2*Cm1;

% Scaled weight
Cw0 = W/(1/2*rho*u0^2*S);

% Dimensionalized stability derivatves
Xu = rho*u0*S*Cw0*sin(thta0) + 1/2*rho*u0*S*Cxu;
Xw = 1/2*rho*u0*S*Cxalpha;
Xq = 1/4*rho*u0*c*S*Cxq;
Xwdot = 1/4*rho*c*S*Cxalphadot;

Zu = -rho*u0*S*Cw0*cos(thta0) + 1/2*rho*u0*S*Czu;
Zw = 1/2*rho*u0*S*Czalpha;
Zq = 1/4*rho*u0*c*S*Czq;
Zwdot = 1/4*rho*c*S*Czalphadot;

Mu = 1/2*rho*u0*c*S*Cmu;
Mw = 1/2*rho*u0*c*S*Cmalpha;
Mq = 1/4*rho*u0*c^2*S*Cmq;
Mwdot = 1/4*rho*c^2*S*Cmalphadot;

% Matrix creation
A1 = [Xu/m Xw/m 0 -g*cos(thta0)];
A2 = [Zu/(m-Zwdot) Zw/(m-Zwdot) (Zq+m*u0)/(m-Zwdot) -m*g*sin(thta0)/(m-Zwdot)];
A3 = 1/Iy * [Mu + Mwdot*Zw/(m-Zwdot) Mw + Mwdot*Zw/(m-Zwdot)...
    Mq + Mwdot*(Zq + m*u0)/(m - Zwdot) -Mwdot*m*g*sin(thta0)/(m-Zwdot)];
A4 = [0 0 1 0];
A_long = [A1; A2; A3; A4];
end
