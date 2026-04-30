function [PHI] = continuousFluence(rs, rd, optProp)
% Giles Blaney Ph.D. Spring 2023
% [PHI] = continuousFluence(rs, rd, optProp)
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
%   PHI     - Fluence. (1/mm^2)

    arguments
        rs (:,3) double;
        rd (:,3) double;

        optProp struct = [];
    end

    PHI=complexFluence(rs, rd, 0, optProp);
end