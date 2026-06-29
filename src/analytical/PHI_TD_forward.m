function [PHI] = PHI_TD_forward(mua, musp, r, t, v)
% PHI_TD_forward Time-domain fluence in an infinite medium.
%
% [PHI] = PHI_TD_forward(mua, musp, r, t, v)
%
% Written by Giles Blaney, Ph.D. Winter 2023
%
% Inputs:
%   mua     - Absorption [1/mm]
%   musp    - Reduced scattering [1/mm]
%   r       - Position vector [nr x 3] [mm]
%   t       - Time after impulse [1 x nt] [sec]
%   v       - Speed of light in medium [mm/sec]
%   Note: r and t should have orthogonal dims
% 
% Outputs:
%   PHI     - Time-domain fluence [1/(mm^2 sec)]
    
    if nargin<=3
        t=(linspace(0, 5e-9, 1000)); %sec
    end
    if nargin<=4
        c=2.99792458e11; %mm/sec
        v=c/1.4;
    end

    rho=sqrt(r(:, 1).^2+r(:, 2).^2);
    z=r(:, 3);
    z0=-1/musp;
    
    PHI=v*((3*(musp+mua))./(4*pi*v*t)).^(3/2).*...
        (exp(-(3*(musp+mua)*((z+z0).^2+rho.^2))./(4*v*t))-...
        exp(-(3*(musp+mua)*((z-z0).^2+rho.^2))./(4*v*t))).*...
        exp(-mua*v*t);
    
    PHI(:, t<=0)=0;
    
end

