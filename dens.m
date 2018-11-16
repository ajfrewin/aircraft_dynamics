function rho = dens(h)
% Calculates atmospheric density under standard atmosphere model
% metric units
if h<=11000
    rho = 1.225 * (temp(h)/288.16)^(-1 - 9.8 / (-6.5e-3 * 287));
else
    rho = dens(11000) * exp(-9.8 * (h - 11000) / (287 * temp(11000)));
end
end
