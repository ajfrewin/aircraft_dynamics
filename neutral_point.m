function h_np = neutral_point(a, at, S, St, l, c, de_dalpha)

% calculates neutral point of aircraft given lift curve slopes and 
% geometric properties

h_ac = .25; % assume aerodynamic center is quater-chord
Vh = l/c * St/S; % tail volume


a_bar = a + at*(1 - de_dalpha)*St/S;
h_np = h_ac + at/a_bar * (1 - de_dalpha) * Vh;

end
