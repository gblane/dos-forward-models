function en = zeroOrdBesselRoots(nRoots)
% zeroOrdBesselRoots Return roots of the zeroth-order Bessel function J0.
%
% en = zeroOrdBesselRoots(nRoots)
%
% Written by Giles Blaney, Ph.D. maintenance update (2026)
%
% Inputs:
%   nRoots - Number of positive J0 roots to return. Default: 2000.
%
% Outputs:
%   en     - Column vector containing the first nRoots positive roots.
%
% Notes:
%   This replaces the historical dependency on zeroOrdBesselRoots.mat, which
%   is not present in the maintained sibling repositories.

    if nargin < 1 || isempty(nRoots)
        nRoots = 2000;
    end
    if nRoots < 1 || nRoots ~= fix(nRoots)
        error('zeroOrdBesselRoots:InvalidNRoots', ...
            'nRoots must be a positive integer.');
    end

    persistent cachedRoots
    if numel(cachedRoots) < nRoots
        startInd = numel(cachedRoots) + 1;
        cachedRoots(nRoots, 1) = NaN;
        for i = startInd:nRoots
            cachedRoots(i, 1) = fzero(@(x) besselj(0, x), ...
                [(i - 0.5)*pi, i*pi]);
        end
    end

    en = cachedRoots(1:nRoots);
end
