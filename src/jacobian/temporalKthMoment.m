function [tkMom] = temporalKthMoment(rs, rd, k, optProp)
% Giles Blaney Ph.D. Spring 2023
% [R] = temporalReflectance(rs, rd, tns, optProp)
% Inputs:
%   rs      - Source corrdinates. (mm)
%   rd      - Detector corrdinates. (mm)
%   k       - (OPTIONAL; default=1)Moment order.
%   optProp - (OPTIONAL) Struct of optical properties with the following
%             fields:
%                nin  - (default=1.4) Index of refraction inside. (-)
%                nout - (default=1) Index of refraction outside. (-)
%                musp - (default=1.2 1/mm) Reduced scattering. (1/mm)
%                mua  - (default=0.01 1/mm) Absorption. (1/mm)
% Outputs:
%   tkMom     - kth Momment of t; <t^k>. (ps^k)
    
    arguments
        rs (:,3) double; %mm
        rd (:,3) double; %mm

        k (1,1) double = 1; 

        optProp struct = [];
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
    
    if size(rs, 1)>1 && size(rd, 1)>1
        error('Can not use multiple sources and multiple detectors');
    end

    nin=optProp.nin;
    nout=optProp.nout;
    musp=optProp.musp; %1/mm
    mua=optProp.mua; %1/mm
    
    x0=rs(:, 1); %mm
    y0=rs(:, 2); %mm
    z0=rs(:, 3); %mm

    c=0.299792458; %mm/ps
    v=c/nin;

    A=n2A(nin, nout);
    D=1/(3*musp); %mm
    zb=-2*A*D; %mm

    mueff=sqrt(mua/D);

    rsp=[x0, y0, -z0+2*zb]; %mm

    r1=vecnorm(rd-rs, 2, 2);
    r2=vecnorm(rd-rsp, 2, 2);
    
    R_C=continuousReflectance(rs, rd, optProp);

    switch k
        case 1
            tkMom=...
                ((z0./r1).*exp(-mueff*r1)+((z0-2*zb)./r2).*exp(-mueff*r2))./...
                (8*pi*v*D*R_C);
        case 2
            tkMom=...
                (z0.*exp(-mueff*r1)+(z0-2*zb).*exp(-mueff*r2))./...
                (16*pi*(v*D)^2*mueff*R_C);
        case 3
            tkMom=...
                (z0.*(r1+1/mueff).*exp(-mueff*r1)+...
                (z0-2*zb).*(r2+1/mueff).*exp(-mueff*r2))./...
                (32*pi*D^2*v^3*mua*R_C);
        case 4
            tkMom=...
                (z0.*(r1.^2+3*r1/mueff+3/mueff^2).*exp(-mueff*r1)+...
                (z0-2*zb).*(r2.^2+3*r2/mueff+3/mueff^2).*exp(-mueff*r2))./...
                (64*pi*D^(5/2)*mua^(3/2)*v^4*R_C);
        otherwise
            error('k>4 not supported');
    end
end