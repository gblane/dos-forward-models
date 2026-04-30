function [L, R] = continuousTotPathLen(rs, rd, optProp)
% Giles Blaney Ph.D. Spring 2023
% [L, R] = continuousTotPathLen(rs, rd, optProp)
% Inputs:
%   rs      - Source corrdinates. (mm)
%   rd      - Detector corrdinates. (mm)
%   optProp - (OPTIONAL) Struct of optical properties with the following
%             fields:
%                nin  - Index of refraction inside. (-)
%                nout - Index of refraction outside. (-)
%                musp - Reduced scattering. (1/mm)
%                mua  - Absorption. (1/mm)
% Outputs:
%   L       - Total pathlength. (mm)
%   R       - Reflectance. (1/mm^2)

    arguments
        rs (:,3) double;
        rd (:,3) double;

        optProp struct = [];
    end

    [L, R]=complexTotPathLen(rs, rd, 0, optProp);
end