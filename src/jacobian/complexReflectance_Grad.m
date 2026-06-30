function [R_Grad] = complexReflectance_Grad(rs, rd, omega, optProp)
% complexReflectance_Grad Calculate the gradient of the complex reflectance.
%
% [R_Grad] = complexReflectance_Grad(rs, rd, omega, optProp)
%
% Written by Giles Blaney (Winter 2020; Ph.D. awarded May 2022)
%
% Inputs:
%   rs      - Source coordinates [mm]
%   rd      - Detector coordinates [mm]
%   omega   - Angular modulation frequency [rad/sec]
%   optProp - Struct of optical properties [struct]
%
% Outputs:
%   R_Grad - Complex reflectance gradient [1/mm^3]

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
    
    r1=vecnorm(rs-rd, 2, 2);
    r2=vecnorm(rsp-rd, 2, 2);
    
    r1_Grad=(rs-rd)./r1;
    r2_Grad=(rsp-rd).*[1, 1, -1]./r2;
    
    term1=(((mueff+1./r1).*exp(-mueff.*r1)./r1.^2+...
        (mueff+1./r2).*exp(-mueff.*r2)./r2.^2)./...
        (4*pi)).*...
        [0, 0, 1];
    
    term2=-((exp(-mueff.*r1)./r1.^4).*...
        z0.*...
        (1+(1+mueff.*r1).*(2+mueff.*r1)).*...
        r1_Grad)./...
        (4*pi);
    
    term3=-((exp(-mueff.*r2)./r2.^4).*...
        (z0-2*zb).*...
        (1+(1+mueff.*r2).*(2+mueff.*r2)).*...
        r2_Grad)./...
        (4*pi);
    
    R_Grad=term1+term2+term3;
end