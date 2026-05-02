function [l_Scat] = complexPartPathLen_Scat(rs, r, rd, V, omega, optProp)
% Giles Blaney Winter 2020
% [l] = complexPartPathLen(rs, rd, omega, optProp)
% Inputs:
%   rs      - Source coordinates. (mm)
%   r       - Center corrdinate of volume. (mm)
%   rd      - Detector coordinates. (mm)
%   V       - Volume. (mm^3)
%   omega   - (OPTIONAL, default=2*pi*1.40625e8 rad/sec) Angular modulation
%             frequency. (rad/sec)
%   optProp - (OPTIONAL) Struct of optical properties with the following
%             fields:
%                nin  - (default=1.4) Index of refraction inside. (-)
%                nout - (default=1) Index of refraction outside. (-)
%                musp - (default=1.2 1/mm) Reduced scattering. (1/mm)
%                mua  - (default=0.01 1/mm) Absorption. (1/mm)
% Outputs:
%   l_Scat  - Complex partial pathlength for scattering. (mm)

    if nargin<=4
        fmod=1.40625e8; %Hz
        omega=2*pi*fmod; %rad/sec
        
        optProp.nin=1.4;
        optProp.nout=1;
        optProp.musp=1.2; %1/mm
        optProp.mua=0.01; %1/mm
    end
    
    PHIrs_r_Grad=complexFluence_Grad(rs, r, omega, optProp);
    Rr_rd_Grad=complexReflectance_Grad(r, rd, omega, optProp);
    Rrs_rd=complexReflectance(rs, rd, omega, optProp);
    
    dotProd=dot(conj(PHIrs_r_Grad), Rr_rd_Grad, 2);
    
    l_Scat=-(dotProd.*V)./(3.*optProp.musp.^2.*Rrs_rd);
end