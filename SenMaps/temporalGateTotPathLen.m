function [L] = temporalGateTotPathLen(rs, rd, tg, optProp, NVA)
% Giles Blaney Ph.D. Spring 2023
% [L] = temporalGateTotPathLen(rs, rd, tg, optProp)
% Inputs:
%   rs      - Source corrdinates. (mm)
%   rd      - Detector corrdinates. (mm)
%   tg      - Gate start and end time. (ps)
%   optProp - (OPTIONAL) Struct of optical properties with the following
%             fields:
%                nin  - (default=1.4) Index of refraction inside. (-)
%                nout - (default=1) Index of refraction outside. (-)
%                musp - (default=1.2 1/mm) Reduced scattering. (1/mm)
%                mua  - (default=0.01 1/mm) Absorption. (1/mm)
%   Name Value Arguments:
%           - 'conv_t' (default=10e3 ps): Time window for convolution. (ps)
%           - 'conv_dt' (default=1 ps): Time step for convolution. (ps)
%           - 'simTyp' (default='DT'): String to switch between 'DT' and
%               'MC' simulation type
% Outputs:
%   L       - Total pathlength of gated t. (mm)
    
    arguments
        rs (:,3) double; %mm
        rd (:,3) double; %mm
        tg (2,:) double; %ps

        optProp struct = [];

        NVA.conv_t (1,1) double = 10e3; %ps
        NVA.conv_dt (1,1) = 1; %ps

        NVA.simTyp (1,:) string = 'DT';
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

    c=0.299792458; %mm/ps
    v=c/optProp.nin;

    t=-NVA.conv_t:NVA.conv_dt:NVA.conv_t; %ps
    
    Rsd_t=temporalReflectance(rs, rd, t, optProp, 'simTyp', NVA.simTyp);

    L=NaN(size(Rsd_t, 1), size(tg, 2));
    for j=1:size(tg, 2)
        [~, i1]=min(abs(t-tg(1, j)));
        [~, i2]=min(abs(t-tg(2, j)));
        
        num=trapz(t(i1:i2), v*t(i1:i2).*Rsd_t(:, i1:i2), 2);
        Rsd_g=trapz(t(i1:i2), Rsd_t(:, i1:i2), 2);

        L(:, j)=num./Rsd_g;
    end
end