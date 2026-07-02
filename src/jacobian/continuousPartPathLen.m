function [l] = continuousPartPathLen(rs, r, rd, V, optProp)
% continuousPartPathLen Calculate continuous-wave (CW) partial pathlength.
%
% [l] = continuousPartPathLen(rs, r, rd, V, optProp)
%
% Written by Giles Blaney, Ph.D. Spring 2023
%
% Inputs:
%   rs      - Source coordinates [mm]
%   r       - Center coordinate of volume [mm]
%   rd      - Detector coordinates [mm]
%   V       - Volume [mm^3]
%   optProp - Struct of optical properties [struct]
%
% Outputs:
%   l - Partial pathlength [mm]

    arguments
        rs (:,3) double;
        r (:,3) double;
        rd (:,3) double;
        V (1,1) double;

        optProp struct = [];
    end

    l = complexPartPathLen(rs, r, rd, V, 0, optProp);
end