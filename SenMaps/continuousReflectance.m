function [R] = continuousReflectance(rs, rd, optProp)
% Giles Blaney Ph.D. Spring 2023
% [R] = continuousReflectance(rs, rd, optProp)
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
%   R       - Reflectance. (1/mm^2)

    arguments
        rs (:,3) double;
        rd (:,3) double;

        optProp struct = [];
    end

    R=complexReflectance(rs, rd, 0, optProp);
end