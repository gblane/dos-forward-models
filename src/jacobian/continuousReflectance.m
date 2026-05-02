function [R] = continuousReflectance(rs, rd, optProp)
% continuousReflectance Calculate continuous-wave (CW) reflectance.
%
% [R] = continuousReflectance(rs, rd, optProp)
%
% Written by Giles Blaney, Ph.D. Spring 2023
%
% Inputs:
%   rs      - Source coordinates [mm]
%   rd      - Detector coordinates [mm]
%   optProp - Struct of optical properties [struct]
%
% Outputs:
%   R - Reflectance [1/mm^2]

    arguments
        rs (:,3) double;
        rd (:,3) double;

        optProp struct = [];
    end

    R=complexReflectance(rs, rd, 0, optProp);
end