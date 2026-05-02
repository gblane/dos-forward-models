function [L_Scat, R] = complexTotPathLen_Scat(rs, rd, omega, optProp)
% Giles Blaney Winter 2020
% [L_Scat, R] = complexTotPathLen_Scat(rs, rd, omega, optProp)
% Inputs:
%   rs      - Source corrdinates. (mm)
%   rd      - Detector corrdinates. (mm)
%   omega   - (OPTIONAL, default=2*pi*1.40625e8 rad/sec) Angular modulation
%             frequecy. (rad/sec)
%   optProp - (OPTIONAL) Struct of optical properties with the following
%             fields:
%                nin  - (default=1.4) Index of refraction inside. (-)
%                nout - (default=1) Index of refraction outside. (-)
%                musp - (default=1.2 1/mm) Reduced scattering. (1/mm)
%                mua  - (default=0.01 1/mm) Absorption. (1/mm)
% Outputs:
%   L_Scat  - Complex total pathlength for scattering. (mm)
%   R       - Complex reflectance. (1/mm^2)

    if nargin<=2
        fmod=1.40625e8; %Hz
        omega=2*pi*fmod; %rad/sec
        
        optProp.nin=1.4;
        optProp.nout=1;
        optProp.musp=1.2; %1/mm
        optProp.mua=0.01; %1/mm
        
        warning(['Default optical properties used, this may be inconsistent'...
            ' with the musp used for source depth']);
    end

    if size(rs, 1)>1 && size(rd, 1)>1
        error('Can not use multiple sources and multiple detectors');
    end
    
    dmusp_numDer=1e-12*optProp.musp;
    optProp1=optProp;
    optProp1.musp=optProp.musp+dmusp_numDer;
    
    R=complexReflectance(rs, rd, omega, optProp);
    R1=complexReflectance(rs, rd, omega, optProp1);
    
    L_Scat=-(1./R).*((R1-R)./dmusp_numDer);
    
end