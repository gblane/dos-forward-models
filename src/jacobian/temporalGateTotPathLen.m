function [L] = temporalGateTotPathLen(rs, rd, tg, optProp, NVA)
% temporalGateTotPathLen Calculate total pathlength for a temporal gate.
%
% [L] = temporalGateTotPathLen(rs, rd, tg, optProp, NVA)
%
% Written by Giles Blaney, Ph.D. Spring 2023
%
% Inputs:
%   rs      - Source coordinates [mm]
%   rd      - Detector coordinates [mm]
%   tg      - Gate start and end time [ps]
%   optProp - Struct of optical properties [struct]
%   NVA     - Name-Value Arguments [struct]
%
% Outputs:
%   L - Total pathlength of gated time [mm]
    
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