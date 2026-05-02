function [V] = temporalVar(rs, rd, optProp)
% Giles Blaney Ph.D. Spring 2023
% [V] = temporalVar(rs, rd, optProp)
% Inputs:
%   rs      - Source corrdinates. (mm)
%   rd      - Detector corrdinates. (mm)
%   optProp - (OPTIONAL) Struct of optical properties with the following
%             fields:
%                nin  - (default=1.4) Index of refraction inside. (-)
%                nout - (default=1) Index of refraction outside. (-)
%                musp - (default=1.2 1/mm) Reduced scattering. (1/mm)
%                mua  - (default=0.01 1/mm) Absorption. (1/mm)
% Outputs:
%   V       - Varinace; <t^2>-<t>^2. (ps^2)
    
    arguments
        rs (:,3) double; %mm
        rd (:,3) double; %mm

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
    
    t1Mom=temporalKthMoment(rs, rd, 1, optProp);
    t2Mom=temporalKthMoment(rs, rd, 2, optProp);

    V=t2Mom-t1Mom.^2;
end