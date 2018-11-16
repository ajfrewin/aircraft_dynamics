function rho = dens_imp(h)
% Calculates atmospheric density under standard atmosphere model
% imperial units
if h<=36089
    rho = 0.0023769 * (temp_ft(h)/288.16)^(-1 - 32.15 / (-1.9812e-3 * 3.0892e3));
else
    rho = dens_imp(36089) * exp(-32.15 * (h - 36089) / (3.0892e3 * temp(36089)));
end
end
