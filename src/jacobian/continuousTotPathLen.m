function [L, R] = continuousTotPathLen(rs, rd, optProp)
% continuousTotPathLen Calculate continuous-wave (CW) total pathlength.
%
% [L, R] = continuousTotPathLen(rs, rd, optProp)
%
% Written by Giles Blaney, Ph.D. Spring 2023
%
% Inputs:
%   rs      - Source coordinates [mm]
%   rd      - Detector coordinates [mm]
%   optProp - Struct of optical properties [struct]
%
% Outputs:
%   L - Total pathlength [mm]
%   R - Reflectance [1/mm^2]

    arguments
        rs (:,3) double;
        rd (:,3) double;

        optProp struct = [];
    end

    [L, R] = complexTotPathLen(rs, rd, 0, optProp);
end