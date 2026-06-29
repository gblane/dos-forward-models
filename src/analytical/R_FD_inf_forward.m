function [R] = R_FD_inf_forward(rho, mua, musp, omega, v)
% R_FD_inf_forward Frequency-domain reflectance from an infinite medium.
%
% [R] = R_FD_inf_forward(rho, mua, musp, omega, v)
%
% Written by Giles Blaney, Ph.D. Spring 2019
%
% Inputs:
%   rho     - Source detector distance [mm]
%   mua     - Absorption [1/mm]
%   musp    - Reduced scattering [1/mm]
%   omega   - Angular modulation frequency [rad/sec]
%   v       - Speed of light in medium [mm/sec]
%                
% Outputs:
%   R       - Complex reflectance [1/mm^2]

% Equation 12.2 in Bigio and Fantini
    
    R=((3*(musp+mua))/(4*pi))*...
        exp(-rho*sqrt(3*(mua+musp)*(mua-1i*omega/v)))./rho;

end

