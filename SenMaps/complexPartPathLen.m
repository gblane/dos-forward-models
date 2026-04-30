function [l] = complexPartPathLen(rs, r, rd, V, omega, optProp)
% Giles Blaney Spring 2019
% Updated Giles Blaney Ph.D. Spring 2023
% [l] = complexPartPathLen(rs, r, rd, V, omega, optProp)
% Inputs:
%   rs      - Source corrdinates. (mm)
%   r       - Center corrdinate of volume. (mm)
%   rd      - Detector corrdinates. (mm)
%   V       - Volume. (mm^3)
%   omega   - (OPTIONAL, default=2*pi*1.40625e8 rad/sec) Angular modulation
%             frequecy. (rad/sec)
%   optProp - (OPTIONAL) Struct of optical properties with the following
%             fields:
%                nin  - (default=1.4) Index of refraction inside. (-)
%                nout - (default=1) Index of refraction outside. (-)
%                musp - (default=1.2 1/mm) Reduced scattering. (1/mm)
%                mua  - (default=0.01 1/mm) Absorption. (1/mm)
% Outputs:
%   l       - Complex partial pathlength. (mm)

    arguments
        rs (:,3) double;
        r (:,3) double;
        rd (:,3) double;
        V (1,1) double;
        
        omega (1,1) = [];
        optProp struct = [];
    end
    
    if isempty(omega)
        fmod=1.40625e8; %Hz
        omega=2*pi*fmod; %rad/sec
    end
        
    if isempty(optProp)
        clear optProp;
        
        optProp.nin=1.4;
        optProp.nout=1;
        optProp.musp=1.2; %1/mm
        optProp.mua=0.01; %1/mm
        
        warning(['Default optical properties used, this may be inconsistent'...
            ' with the musp used for source depth']);
    end
    
    PHIrs_r=complexFluence(rs, r, omega, optProp);
    Rr_rd=complexReflectance(r, rd, omega, optProp);
    Rrs_rd=complexReflectance(rs, rd, omega, optProp);
    
    l=(PHIrs_r.*Rr_rd.*V)./Rrs_rd;
end