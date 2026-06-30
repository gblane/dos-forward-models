function P0_LED = LEDspec_func(lam, lamPK, lamFWHM)
% P0_LED = LEDspec_func(lam, lamPK, lamFWHM)
%   Normalized Gaussian LED emission spectrum (peak-normalized to one). Models a
%   light-emitting-diode source as a Gaussian centered at lamPK with full width
%   at half maximum lamFWHM, sampled on the wavelength grid lam.
%
%   Inputs:
%       lam     - Wavelength grid [nm]
%       lamPK   - Peak (center) wavelength [nm]
%       lamFWHM - Full width at half maximum [nm]. If lamFWHM <= 0 the spectrum
%                 collapses to a unit impulse at the grid point nearest lamPK
%                 (the monochromatic / laser-diode limit).
%   Output:
%       P0_LED  - Emission spectrum sampled at lam, with max(P0_LED) == 1 [-]

    if lamFWHM <= 0
        % Monochromatic / laser-diode limit: unit impulse at the nearest grid point
        P0_LED = zeros(size(lam));
        [~, iPk] = min(abs(lam - lamPK));
        P0_LED(iPk) = 1;
        return
    end

    % Gaussian (FWHM -> standard deviation), normalized to a unit peak
    P0_LED = normpdf(lam, lamPK, lamFWHM/(2*sqrt(2*log(2))));
    P0_LED = P0_LED/max(P0_LED);

end
