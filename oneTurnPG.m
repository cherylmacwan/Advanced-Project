function PG = oneTurnPG(wavelength, ht, hr, x, y, D, exact, polarization, eps)
% Simplified model of PG after turning a corner 
% ht = transmitter Height [m]
% hr = reciever Height [m]
% x = distance to corner [m]
% y = distance from corner to reciever [m]
% D = Diffraction coefficient

% PG = (wavelength./ (4.*pi.*(x)) ).^2 .* (2.*
% sin(2.*pi.*ht.*hr./(wavelength.*(x+y)))).^2 *(4.*D^2).*(x+y)./(x*y);

if(exact)
    PG = exact2RayModel(ht, hr, x+y, polarization, eps, wavelength) *(D^2).*(x+y)./(x.*(y));
else
    PG = PG2R2(wavelength, ht, hr, (x+y))*(D^2).*(x+y)./(x.*(y));
end