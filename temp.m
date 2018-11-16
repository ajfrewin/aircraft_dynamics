function tau = temp(h)
% Calculates atmospheric temp under standard atmosphere model
% metric units
if h<=11000
    tau = 288.16 - 6.5e-3 * h;
else
    tau = 288.16 - 6.5e-3 * 11000;
end

end
