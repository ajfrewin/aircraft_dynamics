function A_lat = lat_stab(Clbeta, Clp, Clr, CYbeta, CYp, CYr, Cnbeta, Cnp, Cnr,...
    m, rho, S, b, thta0, u0, Ix, Iz, Izx)
% Author: Adam Frewin 2018
% Description: Computes lateral stability matrix of a standard fixed-wing
% aircraft, according to the lateral linear model

g = 32.26; % gravity acceleration, ft/s^2
scale = 1/2*rho*u0*S*b;  % dimensional scale factor
Ixprim = (Ix*Iz - Izx^2)/Iz;  % Scaled moments of inertia
Izprim = (Ix*Iz - Izx^2)/Ix;
Izxprim = Izx/(Ix*Iz - Izx^2);

% Dimensionalized Stability derivatives
Yv = scale/b*CYbeta;
Yp = 1/2*scale*CYp;
Yr = 1/2*scale*CYr;

Lv = scale*Clbeta;
Lp = 1/2*scale*b*Clp;
Lr = 1/2*scale*b*Clr;

Nv = scale*Cnbeta;
Np = 1/2*scale*b*Cnp;
Nr = 1/2*scale*b*Cnr;

% Compilation of the matrix
A1 = [Yv/m Yp/m (Yr/m - u0) g*cos(thta0)];
A2 = [(Lv/Ixprim + Izxprim*Nv) (Lp/Ixprim + Izxprim*Np)...
    (Lr/Ixprim + Izxprim*Nr) 0];
A3 = [(Izxprim*Lv + Nv/Izprim) (Izxprim*Lp + Np/Izprim)...
    (Izxprim*Lr + Nr/Izprim) 0];
A4 = [0 1 tan(thta0) 0];

A_lat = [A1; A2; A3; A4];
end