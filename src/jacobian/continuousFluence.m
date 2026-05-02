function [PHI] = continuousFluence(rs, rd, optProp)
% continuousFluence Calculate continuous-wave (CW) fluence.
%
% [PHI] = continuousFluence(rs, rd, optProp)
%
% Written by Giles Blaney, Ph.D. Spring 2023
%
% Inputs:
%   rs      - Source coordinates [mm]
%   rd      - Detector coordinates [mm]
%   optProp - Struct of optical properties [struct]
%
% Outputs:
%   PHI - Fluence [1/mm^2]

    arguments
        rs (:,3) double;
        rd (:,3) double;

        optProp struct = [];
    end

    PHI=complexFluence(rs, rd, 0, optProp);
end