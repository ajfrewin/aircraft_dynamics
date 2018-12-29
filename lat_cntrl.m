function B_lat = lat_cntrl(Cydeltaa, Cydeltar, Cldeltaa, Cldeltar, CNdeltaa, CNdeltar,...
    rho, u0, S, b, m, Ix, Iz, Izx)
% Author: Adam Frewin 2018
% Description: Computes the lateral control matrix for fixed-wing aircraft
% under lateral-linear model

% Scaled moments of inertia
Ixprim = (Ix*Iz - Izx^2)/Iz;
Izprim = (Ix*Iz - Izx^2)/Ix;
Izxprim = Izx / (Ix*Iz - Izx^2);

sc = 1/2*rho*u0^2*S;  % dimensional scale factor

% Dimensionalized 
Ydeltaa = Cydeltaa*sc;
Ydeltar = Cydeltar*sc;
Ldeltaa = Cldeltaa*sc*b;
Ldeltar = Cldeltar*sc*b;
Ndeltaa = CNdeltaa*sc*b;
Ndeltar = CNdeltar*sc*b;

% Matrix Creation
B1 = [Ydeltaa/m Ydeltar/m];
B2 = [Ldeltaa/Ixprim + Izxprim*Ndeltaa Ldeltar/Ixprim + Izxprim*Ndeltar];
B3 = [Ndeltaa/Izprim + Izxprim*Ldeltaa Ndeltar/Izprim + Izxprim*Ldeltar];
B4 = [0 0];

B_lat = [B1; B2; B3; B4];