function [ PG ] = exact2RayModel(ht, hr, R, polarization, eps, wavelength)
% Function to compute the exact 2 ray model

r1 = ((ht-hr)^2 + R.^2).^0.5;
r2 = ((ht+hr)^2 + R.^2).^0.5;
theta = atan(R/(ht+hr));

if(polarization == 0)
    a = 1/eps;
elseif(polarization == 1)
    a = 1;
else
   error('Error: Invalid Polarization') 
end

gamma = (cos(theta)-a*(eps - sin(theta).^2).^0.5) ...
    ./ (cos(theta)+a*(eps - sin(theta).^2).^0.5);


PG = (wavelength/(4*pi))^2*abs(exp(-1i*2*pi/wavelength*r1)./r1 + gamma.*exp(-1i*2*pi/wavelength*r2)./r2).^2; 


end

