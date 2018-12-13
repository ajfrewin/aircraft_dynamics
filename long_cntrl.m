function B_long = long_cntrl(CDdeltae, CTxu, CTx1, CLdeltae, Cmdeltae, CmTu,...
    CmT1, rho, u0, S, m, c, Zwdot, Mwdot, Iy)

Cxdeltae = -CDdeltae;
Cxdeltap = CTxu + 2*CTx1;
Czdeltae = -CLdeltae;
Czdeltap = 0;
Cmdeltap = CmTu + 2*CmT1;

sc = 1/2*rho*u0^2*S;
Xdeltae = sc*Cxdeltae;
Xdeltap = sc*Cxdeltap;

Zdeltae = Czdeltae*sc;
Zdeltap = Czdeltap*sc;

Mdeltae = Cmdeltae*sc*c;
Mdeltap = Cmdeltap*sc*c;

B1 = [Xdeltae/m Xdeltap/m];
B2 = [Zdeltae/(m-Zwdot) Zdeltap/(m-Zwdot)];
B3 = [Mdeltae/Iy + Mwdot*Zdeltae/(Iy*(m-Zwdot)) Mdeltap/Iy + Mwdot*Zdeltap/(Iy*(m-Zwdot))];
B4 = [0 0];
B_long = [B1; B2; B3; B4];
end
