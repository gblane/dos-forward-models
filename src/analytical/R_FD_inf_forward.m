function [R] = R_FD_inf_forward(rho, mua, musp, omega, v)
% [R] = R_FD_forward(rho, mua, musp, omega)
% Giles Blaney Spring 2019
% Inputs:
%   rs      - Source detector distance. (mm)
%   rd      - Detector coordinates. (mm)
%   musp    - Reduced scattering. (1/mm)
%   mua     - Absorption. (1/mm)
%   omega   - Angular modulation frequency. (rad/sec)
%                
% Outputs:
%   R       - Complex reflectance. (1/mm^2)

% Equation 12.2 in Bigio and Fantini
    
    R=((3*(musp+mua))/(4*pi))*...
        exp(-rho*sqrt(3*(mua+musp)*(mua-1i*omega/v)))./rho;

end

