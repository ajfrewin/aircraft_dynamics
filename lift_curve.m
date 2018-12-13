function a = lift_curve(AR, Mach, sweep)
if AR<4
    k = 1 + AR*(1.87 - 2.33e-4*sweep)/100;
else
    k = 1 + (8.2 - 2.3*sweep - AR*(.22-.153*sweep))/100;
end

a = 2*pi*AR./(2+sqrt(AR^2.*(1-Mach.^2)./k^2 .* (1 + tan(sweep)^2./(1 - Mach.^2))+4));
end
