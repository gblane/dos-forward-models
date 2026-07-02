function [L] = temporalVarTotPathLen(rs, rd, optProp)
% temporalVarTotPathLen Calculate total pathlength of the temporal variance.
%
% [L] = temporalVarTotPathLen(rs, rd, optProp)
%
% Written by Giles Blaney, Ph.D. Spring 2023
%
% Inputs:
%   rs      - Source coordinates [mm]
%   rd      - Detector coordinates [mm]
%   optProp - Struct of optical properties [struct]
%
% Outputs:
%   L - Total pathlength of variance [mm]

    arguments
        rs (:,3) double; % mm
        rd (:,3) double; % mm

        optProp struct = [];
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

    t1Mom = temporalKthMoment(rs, rd, 1, optProp);
    t2Mom = temporalKthMoment(rs, rd, 2, optProp);
    V = t2Mom-t1Mom.^2;

    Lt1 = temporalKthMomTotPathLen(rs, rd, 1, optProp);
    Lt2 = temporalKthMomTotPathLen(rs, rd, 2, optProp);

    L = (t2Mom.*Lt2-2*t1Mom.^2.*Lt1)./V;
end