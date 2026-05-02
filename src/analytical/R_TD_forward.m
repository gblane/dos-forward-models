function [R] = R_TD_forward(rho, mua, musp, t, v, z0)
% [R] = R_FD_forward(rho, mua, musp, t, v)
% Giles Blaney Summer 2021
% Uses Zero-Boundry (Equ. 11.18 in Bigio & Fantini)
% Inputs:
%   rho     - Source detector distance. (mm)
%   musp    - Reduced scattering. (1/mm)
%   mua     - Absorption. (1/mm)
%   t       - Time after impulse. (sec)
%   v       - Speed of light in medium. (mm/sec)
%   z0      - Iso-Source Depth
%   Note: rho and t should have orthogonal dims
% 
% Outputs:
%   R       - Complex reflectance. (1/(mm^2 sec))
    
    if nargin<=3
        t=(linspace(0, 5e-9, 1000)); %sec
    end
    if nargin<=4
        c=2.99792458e11; %mm/sec
        v=c/1.4;
    end
    if nargin<=5
        z0=1/musp;
    end
    
    R=((3*(musp+mua))/(4*pi*v))^(3/2)*...
        (z0./(t.^(5/2))).*...
        exp(-(3*(musp+mua)*(z0.^2+rho.^2))./(4*v*t)-mua*v*t);
    
    R(t<=0)=0;
    
end

