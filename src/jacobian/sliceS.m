function [Splane, plotParams] = sliceS(S, params, slAxis, slValue)
% sliceS Slice a 3D sensitivity array for 2D plotting.
%
% [Splane, plotParams] = sliceS(S, params, slAxis, slValue)
%
% Written by Giles Blaney, Ph.D. Spring 2023
%
% Inputs:
%   S       - 3D sensitivity array [unitless]
%   params  - Parameters structure containing axis vectors [struct]
%   slAxis  - Axis to slice along ('x', 'y', or 'z') [string]
%   slValue - Value along the slice axis [mm]
%
% Outputs:
%   Splane     - 2D sliced sensitivity array [unitless]
%   plotParams - Plotting parameters and axis vectors [struct]

    slAxis = lower(slAxis);
    [~, ind] = min(abs(params.(slAxis)-slValue));

    switch slAxis
        case "x"
            Splane = squeeze(S(ind, :, :)).';
            plotParams.vertAx = params.z;
            plotParams.vertNm = "$z$ (mm)";
            plotParams.horzAx = params.y;
            plotParams.horzNm = "$y$ (mm)";
        case "y"
            Splane = squeeze(S(:, ind, :)).';
            plotParams.vertAx = params.z;
            plotParams.vertNm = "$z$ (mm)";
            plotParams.horzAx = params.x;
            plotParams.horzNm = "$x$ (mm)";
        case "z"
            Splane = S(:, :, ind).';
            plotParams.vertAx = params.y;
            plotParams.vertNm = "$y$ (mm)";
            plotParams.horzAx = params.x;
            plotParams.horzNm = "$x$ (mm)";
        otherwise
            error("Unknown axis %s", slAxis);
    end
end