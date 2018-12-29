function alpha_0=static_stability(h, h_ac, a, at, S, St, l, c, de_dalpha,...
    Cmac, Clow, it, e0, plot_vals)
% Author: Adam Frewin 2018
% Description: determines if aircraft is statically stable and if so
% returns the Angle of attack at which it will achieve steady level flight

fprintf('CG = ' + string(h*c*12) +...
    ' inches\nCL_0,w = ' + string(Clow) + '\n');
Vh = l/c * St/S;
Cm0 = Cmac + Clow*(h-h_ac) - at*(it-e0)*Vh * (1 - (h - h_ac)*c/l);
fprintf('Cm0 = ' + string(Cm0) + '\n');
dCm_dalpha = (a + at * St/S * (1 - de_dalpha))*(h-h_ac) - at*Vh*(1-de_dalpha);
fprintf('Cm_alpha = ' + string(dCm_dalpha) + '\n');

% for plots
alpha_span = -5:.1:10;
alpha_span = alpha_span*pi/180;

if(plot_vals)
    plot(alpha_span*180/pi, Cm0 + dCm_dalpha*alpha_span);
end
if dCm_dalpha>=0
    fprintf('Aircraft is not statically stable\n');
    alpha_0 = NaN;
else
    alpha_0 = -Cm0 / dCm_dalpha;
    fprintf('Trim alpha = ' + string(alpha_0*180/pi) + ' degrees\n \n');
end
end
