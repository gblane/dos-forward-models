function [R] = R_TD_forward(rho, mua, musp, t, v, z0)
% R_TD_forward Time-domain reflectance from a semi-infinite medium.
%
% [R] = R_TD_forward(rho, mua, musp, t, v, z0)
%
% Written by Giles Blaney (Summer 2021; Ph.D. awarded May 2022)
%
% Inputs:
%   rho     - Source detector distance [mm]
%   mua     - Absorption [1/mm]
%   musp    - Reduced scattering [1/mm]
%   t       - Time after impulse [sec]
%   v       - Speed of light in medium [mm/sec]
%   z0      - Iso-Source Depth [mm]
%   Note: rho and t should have orthogonal dims
% 
% Outputs:
%   R       - Time-domain reflectance [1/(mm^2 sec)]
    
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

