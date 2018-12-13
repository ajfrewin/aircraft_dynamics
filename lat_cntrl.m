function B_lat = lat_cntrl(Cydeltaa, Cydeltar, Cldeltaa, Cldeltar, CNdeltaa, CNdeltar,...
    rho, u0, S, b, m, Ix, Iz, Izx)

Ixprim = (Ix*Iz - Izx^2)/Iz;
Izprim = (Ix*Iz - Izx^2)/Ix;
Izxprim = Izx / (Ix*Iz - Izx^2);
sc = 1/2*rho*u0^2*S;

Ydeltaa = Cydeltaa*sc;
Ydeltar = Cydeltar*sc;
Ldeltaa = Cldeltaa*sc*b;
Ldeltar = Cldeltar*sc*b;
Ndeltaa = CNdeltaa*sc*b;
Ndeltar = CNdeltar*sc*b;

B1 = [Ydeltaa/m Ydeltar/m];
B2 = [Ldeltaa/Ixprim + Izxprim*Ndeltaa Ldeltar/Ixprim + Izxprim*Ndeltar];
B3 = [Ndeltaa/Izprim + Izxprim*Ldeltaa Ndeltar/Izprim + Izxprim*Ldeltar];
B4 = [0 0];

B_lat = [B1; B2; B3; B4];