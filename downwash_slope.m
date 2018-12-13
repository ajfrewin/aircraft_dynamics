function de_dalpha = downwash_slope(AR, Mach, ltv, taper, span, sweep)
de_dalpha = 4.44*sqrt(1-Mach.^2).*((1/AR - 1/(1 + AR^1.7))*(10-3*taper)/7.*...
    (1 - ltv/span)/((2*lt/b)^.33)*sqrt(cos(sweep))).^1.19;
end
