function [l] = temporalVarPartPathLen(rs, r, rd, Vol, optProp, NVA)
% temporalVarPartPathLen Calculate partial pathlength of the temporal variance.
%
% [l] = temporalVarPartPathLen(rs, r, rd, Vol, optProp, NVA)
%
% Written by Giles Blaney, Ph.D. Spring 2023
%
% Inputs:
%   rs      - Source coordinates [mm]
%   r       - Center coordinate of volume [mm]
%   rd      - Detector coordinates [mm]
%   Vol     - Volume [mm^3]
%   optProp - Struct of optical properties [struct]
%   NVA     - Name-Value Arguments [struct]
%
% Outputs:
%   l - Partial pathlength of variance [mm]
%
% Shared-repo dependencies:
%   struct2pairs is provided by ../my-matlab.

    arguments
        rs (1,3) double; % mm
        r (:,3) double; % mm
        rd (1,3) double; % mm
        Vol (1,1) double; % mm^3

        optProp struct = [];

        NVA.conv_t (1,1) double = 10e3; % ps
        NVA.conv_dt (1,1) = 1; % ps
        NVA.usePar (1,1) logical = true;
        NVA.FFTconv (1,1) logical = true;
    end

    if isempty(optProp)
        clear optProp;

        optProp.nin = 1.4;
        optProp.nout = 1;
        optProp.musp = 1.2; % 1/mm
        optProp.mua = 0.01; % 1/mm

        warning(['Default optical properties used, this may be inconsistent'...
            ' with the musp used for source depth']);
    end

    if size(rs, 1) > 1 && size(rd, 1) > 1
        error("Can not use multiple sources and multiple detectors");
    end

    NVAstruct = struct2pairs(NVA);

    t1Mom = temporalKthMoment(rs, rd, 1, optProp);
    t2Mom = temporalKthMoment(rs, rd, 2, optProp);
    V = t2Mom-t1Mom.^2;

    lt1 = temporalKthMomPartPathLen(rs, r, rd, Vol, 1, optProp, NVAstruct{:});
    lt2 = temporalKthMomPartPathLen(rs, r, rd, Vol, 2, optProp, NVAstruct{:});

    l = (t2Mom.*lt2-2*t1Mom.^2.*lt1)./V;

    l(isnan(l)) = 0;
end
