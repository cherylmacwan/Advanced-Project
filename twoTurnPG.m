function PG = twoTurnPG(wavelength, ht, hr, x, y, z, D, exact, polarization, eps)
% Simplified model of PG after turning a corner 
% ht = transmitter Height [m]
% hr = reciever Height [m]
% x = distance to first intersection from transmitter [m]
% y = distance to second intersection from first intersection [m]
% z = distance from second intersection to reciever [m]
% D = Diffraction coefficient

if(exact)
    PG = exact2RayModel(ht, hr, x+y+z, polarization, eps, wavelength) *(D^4).*(x+y+z)./(x*y*z);
else
    PG = PG2R2(wavelength, ht, hr, (x+y+z))*(D^4).*(x+y+z)./(x*y*z);
end