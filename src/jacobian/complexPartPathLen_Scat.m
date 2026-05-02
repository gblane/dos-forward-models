function [l_Scat] = complexPartPathLen_Scat(rs, r, rd, V, omega, optProp)
% complexPartPathLen_Scat Calculate complex partial pathlength for scattering.
%
% [l_Scat] = complexPartPathLen_Scat(rs, r, rd, V, omega, optProp)
%
% Written by Giles Blaney, Ph.D. Winter 2020
%
% Inputs:
%   rs      - Source coordinates [mm]
%   r       - Center coordinate of volume [mm]
%   rd      - Detector coordinates [mm]
%   V       - Volume [mm^3]
%   omega   - Angular modulation frequency [rad/sec]
%   optProp - Struct of optical properties [struct]
%
% Outputs:
%   l_Scat - Complex partial pathlength for scattering [mm]

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