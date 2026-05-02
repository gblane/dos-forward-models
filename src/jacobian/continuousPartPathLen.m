function [l] = continuousPartPathLen(rs, r, rd, V, optProp)
% Giles Blaney Ph.D. Spring 2023
% [l] = continuousPartPathLen(rs, r, rd, V, optProp)
% Inputs:
%   rs      - Source corrdinates. (mm)
%   r       - Center corrdinate of volume. (mm)
%   rd      - Detector corrdinates. (mm)
%   V       - Volume. (mm^3)
%   optProp - (OPTIONAL) Struct of optical properties with the following
%             fields:
%                nin  - Index of refraction inside. (-)
%                nout - Index of refraction outside. (-)
%                musp - Reduced scattering. (1/mm)
%                mua  - Absorption. (1/mm)
% Outputs:
%   l       - Partial pathlength. (mm)

    arguments
        rs (:,3) double;
        r (:,3) double;
        rd (:,3) double;
        V (1,1) double;

        optProp struct = [];
    end
    
    l=complexPartPathLen(rs, r, rd, V, 0, optProp);
end