function [Splane, plotParams] = sliceS(S, params, slAxis, slValue)
% Giles Blaney Ph.D. Spring 2023
% [Splane, plotParams] = sliceS(S, params, slAxis, slValue)
%
% Slice a 3D S array for 2D plotting

    slAxis=lower(slAxis);
    [~, ind]=min(abs(params.(slAxis)-slValue));

    switch slAxis
        case 'x'
            Splane=squeeze(S(ind, :, :)).';
            plotParams.vertAx=params.z;
            plotParams.vertNm='$z$ (mm)';
            plotParams.horzAx=params.y;
            plotParams.horzNm='$y$ (mm)';
        case 'y'
            Splane=squeeze(S(:, ind, :)).';
            plotParams.vertAx=params.z;
            plotParams.vertNm='$z$ (mm)';
            plotParams.horzAx=params.x;
            plotParams.horzNm='$x$ (mm)';
        case 'z'
            Splane=S(:, :, ind).';
            plotParams.vertAx=params.y;
            plotParams.vertNm='$y$ (mm)';
            plotParams.horzAx=params.x;
            plotParams.horzNm='$x$ (mm)';
        otherwise
            error('Unknown axis %s', slAxis);
    end
end