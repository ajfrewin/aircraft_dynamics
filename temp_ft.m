function tau = temp_ft(h)
% Calculates atmospheric temp under standard atmosphere model
% imperial units
if h<=36089
    tau = 288.16 - 1.9812e-3 * h;
else
    tau = 288.16 - 1.9812e-3 * 11000;
end

end