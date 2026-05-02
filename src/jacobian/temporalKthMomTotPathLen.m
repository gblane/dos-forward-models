function [L] = temporalKthMomTotPathLen(rs, rd, k, optProp)
% Giles Blaney Ph.D. Spring 2023
% [L] = temporalKthMomTotPathLen(rs, rd, k, optProp)
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
%   L       - Total pathlength of kth momment of t. (mm)
    
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

    c=0.299792458; %mm/ps
    v=c/nin;

    t1Mom=temporalKthMoment(rs, rd, 1, optProp); %ps
    tkMom=temporalKthMoment(rs, rd, k, optProp); %ps^k
    tkp1Mom=temporalKthMoment(rs, rd, k+1, optProp); %ps^(k+1)

    L=-(v*(t1Mom.*tkMom-tkp1Mom))./tkMom;
end