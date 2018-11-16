function [alpha_0, delta_e]=static_stability_elevator(h, hac, a, at, ae,...
    S, St, l, c, de_dalpha, Cmac, Clow, it, e0, rho, V, W)

% Determines if aircraft is statically stable and if so
% returns the trim angle of attack

Vh = l/c * St/S; % tail volume

% discretize into Control and Stability Derivatives
% Cl = Cl_0 + Cl_alpha * alpha + Cl_delta * delta_e
Cl_0 = Clow + St/S*at*(it-e0);
fprintf('Cl0 = ' + string(Cl_0) + '\n');

Cl_alpha = a + at*(1 - de_dalpha)*St/S;
fprintf('Cl_alpha = ' + string(Cl_alpha) + '\n');

Cl_delta = St/S * ae;
fprintf('Cl_delta = ' + string(Cl_delta) + '\n');

% Cm = Cm_0 + Cm_alpha * alpha + Cm_delta * delta_e
Cm_0 = Cmac + Clow*(h-hac) - at*(it-e0)*Vh * (1 - (h - hac)*c/l);
fprintf('Cm_0 = ' + string(Cm_0) + '\n');

Cm_alpha = (a + at * St/S * (1 - de_dalpha))*(h-hac) - at*Vh*(1-de_dalpha);
fprintf('Cm_alpha = ' + string(Cm_alpha) + '\n');

Cm_delta = Cl_delta *(h-hac) - ae*Vh;
fprintf('Cm_delta = ' + string(Cm_delta) + '\n');

% set up matrix 
stability_matrix = [Cl_alpha Cl_delta; Cm_alpha Cm_delta];
RHS = [W/(.5*rho*V^2*S) - Cl_0; -Cm_0];
fprintf('Stability matrix : \n');
disp(stability_matrix);

if Cm_alpha>=0
    fprintf('Aircraft is not statically stable\n');
    [alpha_0, delta_e]= NaN;
elseif det(stability_matrix)==0
    fprintf('Warning: Stability matrix is singular, cannot compute using matricies.\n');
    fprintf('Running calculation neglecting elevator: \n');
    alpha_0 = -Cm_0 / Cm_alpha;
    delta_e = NaN;   
    fprintf('Trim alpha = ' + string(alpha_0*180/pi) + ' degrees\n \n');
else
    angles = stability_matrix\RHS;
    alpha_0 = angles(1);
    delta_e = angles(2);
    fprintf('Trim alpha = ' +  string(alpha_0*180/pi) + ' degrees\n');
    fprintf('Trim elevator deflection = ' +  string(delta_e*180/pi) + ' degrees\n \n');
end
fprintf('Alternate calculation:\n');
Cl_trim = W/(.5*rho*V^2*S);
a_bar = a*(1 + at/a * St/S * (1- de_dalpha));
alpha_other = (Cm_0 * Cl_delta + Cm_delta*(Cl_trim - Cl_0))/...
    (a_bar * Cm_delta - Cl_delta*Cm_alpha);
delta_other = -(a_bar * Cm_0 + Cm_alpha *(Cl_trim - Cl_0))/...
    (a_bar*Cm_delta - Cl_delta * Cm_alpha);
fprintf('Trim alpha = ' +  string(alpha_other*180/pi) + ' degrees\n');
fprintf('Trim elevator deflection = ' +  string(delta_other*180/pi) + ' degrees\n \n');
end
