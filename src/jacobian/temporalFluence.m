function [PHI] = temporalFluence(rs, rd, t, optProp)
% Giles Blaney Ph.D. Spring 2023
% [PHI] = temporalFluence(rs, rd, tns, optProp)
% Inputs:
%   rs      - Source corrdinates. (mm)
%   rd      - Detector corrdinates. (mm)
%   t       - Time. (ps)
%   optProp - (OPTIONAL) Struct of optical properties with the following
%             fields:
%                nin  - (default=1.4) Index of refraction inside. (-)
%                nout - (default=1) Index of refraction outside. (-)
%                musp - (default=1.2 1/mm) Reduced scattering. (1/mm)
%                mua  - (default=0.01 1/mm) Absorption. (1/mm)
% Outputs:
%   PHI     - Temporal fluence. (1/(ps mm^2))
    
    arguments
        rs (:,3) double; %mm
        rd (:,3) double; %mm
        t (1,:) double; %ps

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

    rsp=[x0, y0, -z0+2*zb]; %mm

    r1=vecnorm(rd-rs, 2, 2);
    r2=vecnorm(rd-rsp, 2, 2);

    posInds=t>0;
    tPos=t(posInds);
    
    PHI=...
        ((v*exp(-mua*v*tPos))./(4*pi*D*v*tPos).^(3/2)).*...
        (exp(-r1.^2./(4*D*v*tPos))-...
        exp(-r2.^2./(4*D*v*tPos))); %1/(ps mm^2)

    PHI=[zeros(size(PHI, 1), sum(~posInds)), PHI];
end