function it=req_incidence(h, hac, a, at, S, St, l, c, de_dalpha,...
    Cmac, Clow, alpha_t, e0)
% Author: Adam Frewin 2018
% Description: given a target trim alpha and all other values, determines
% required incidence angle of tail

Vh = l/c * St/S; % Tail volume

dCm_dalpha = (a + at*St/S*(1-de_dalpha))*(h-hac) - at*Vh*(1-de_dalpha);

it = (Cmac + Clow*(h-hac) + dCm_dalpha*alpha_t)/(at*Vh*(1-(h-hac)*c/l)) + e0;
end