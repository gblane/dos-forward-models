function [PHI_Grad] = complexFluence_Grad(rs, rd, omega, optProp)
% complexFluence_Grad Calculate the gradient of the complex fluence.
%
% [PHI_Grad] = complexFluence_Grad(rs, rd, omega, optProp)
%
% Written by Giles Blaney, Ph.D. Winter 2020
%
% Inputs:
%   rs      - Source coordinates [mm]
%   rd      - Detector coordinates [mm]
%   omega   - Angular modulation frequency [rad/sec]
%   optProp - Struct of optical properties [struct]
%
% Outputs:
%   PHI_Grad - Complex fluence gradient [1/mm^3]

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
    
    x0=rs(:, 1); %mm
    y0=rs(:, 2); %mm
    z0=rs(:, 3); %mm

    c=2.99792458e11; %mm/sec
    v=c/optProp.nin;

    A=n2A(optProp.nin, optProp.nout);
    D=1/(3*optProp.musp); %mm
    zb=-2*A*D; %mm

    mueff=sqrt(optProp.mua/D-1i*omega/(v*D)); %1/mm

    rsp=[x0, y0, -z0+2*zb]; %mm

    r1=vecnorm(rd-rs, 2, 2);
    r2=vecnorm(rd-rsp, 2, 2);
    
    r1_Grad=(rd-rs)./r1;
    r2_Grad=(rd-rsp)./r2;
    
    PHI_Grad=-(...
        (exp(-mueff.*r1)./r1).*(mueff+1./r1).*r1_Grad-...
        (exp(-mueff.*r2)./r2).*(mueff+1./r2).*r2_Grad)/...
        (4*pi*D);
end