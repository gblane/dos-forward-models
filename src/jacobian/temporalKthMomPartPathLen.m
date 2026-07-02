function [l] = temporalKthMomPartPathLen(rs, r, rd, Vol, k, optProp, NVA)
% temporalKthMomPartPathLen Calculate partial pathlength of the k-th temporal moment.
%
% [l] = temporalKthMomPartPathLen(rs, r, rd, Vol, k, optProp, NVA)
%
% Written by Giles Blaney, Ph.D. Spring 2023
%
% Inputs:
%   rs      - Source coordinates [mm]
%   r       - Center coordinate of volume [mm]
%   rd      - Detector coordinates [mm]
%   Vol     - Volume [mm^3]
%   k       - Moment order [unitless]
%   optProp - Struct of optical properties [struct]
%   NVA     - Name-Value Arguments [struct]
%
% Outputs:
%   l - Partial pathlength of k-th moment of t [mm]

    arguments
        rs (1,3) double; % mm
        r (:,3) double; % mm
        rd (1,3) double; % mm
        Vol (1,1) double; % mm^3

        k (1,1) double = 1;

        optProp struct = [];

        NVA.conv_t (1,1) double = 10e3; % ps
        NVA.conv_dt (1,1) = 1; % ps
        NVA.usePar (1,1) logical = true;
        NVA.FFTconv (1,1) logical = true;
    end

    if isempty(optProp)
        clear optProp;

        optProp.nin = 1.4;
        optProp.nout = 1;
        optProp.musp = 1.2; % 1/mm
        optProp.mua = 0.01; % 1/mm

        warning(['Default optical properties used, this may be inconsistent'...
            ' with the musp used for source depth']);
    end

    if NVA.usePar
        pp = gcp;
    end

    t = -NVA.conv_t:NVA.conv_dt:NVA.conv_t; % ps
    posInds = t > 0;

    tkMom = temporalKthMoment(rs, rd, k, optProp);
    RC = continuousReflectance(rs, rd, optProp);

    if NVA.usePar
        FFTconv = NVA.FFTconv;
        dt = NVA.conv_dt;
        l = NaN(size(r, 1), 1);
        parfor i = 1:size(r, 1)
            lC = continuousPartPathLen(rs, r(i, :), rd, Vol, optProp);

            PHIsi = zeros(size(t));
            Rid = zeros(size(t));
            PHIsi(posInds) = temporalFluence(rs, r(i, :), t(posInds), ...
                optProp);
            Rid(posInds) = temporalReflectance(r(i, :), rd, t(posInds), ...
                optProp);

            if FFTconv
                convPHIsiRid = ifft(fft(PHIsi).*fft(Rid));
                convPHIsiRid((sum(posInds)+1):end) = [];

                l(i) = -(lC*tkMom-(Vol/RC)*...
                    sum(t(posInds).^k.*...
                    (convPHIsiRid*dt),...
                    2))*dt/tkMom;
            else
                l(i) = -(lC*tkMom-(Vol/RC)*...
                    trapz(t, t.^k.*...
                    (conv(PHIsi, Rid, "same")*dt),...
                    2))/tkMom;
            end
        end
    else
        PHIsi = zeros(size(r, 1), length(t));
        Rid = zeros(size(r, 1), length(t));
        lC = continuousPartPathLen(rs, r, rd, Vol, optProp);
        PHIsi(:, posInds) = temporalFluence(rs, r, t(posInds), optProp);
        Rid(:, posInds) = temporalReflectance(r, rd, t(posInds), optProp);
        convPR = NaN(size(r, 1), length(t));
        convPHIsiRid = NaN(size(r, 1), sum(posInds));

        if NVA.FFTconv
            for i = 1:size(r, 1)
                tmp = ifft(fft(PHIsi(i, :)).*fft(Rid(i, :)));
                convPHIsiRid(i, :) = tmp(1:sum(posInds));
            end

            dkMomIdmua = lC*tkMom-(Vol/RC)*...
                sum(...
                t(posInds).^k.*(convPHIsiRid*NVA.conv_dt)*NVA.conv_dt, 2);
        else
            for i = 1:size(r, 1)
                convPR(i, :) = conv( ...
                    PHIsi(i, :), Rid(i, :), "same")*NVA.conv_dt;
            end

            dkMomIdmua = lC*tkMom-(Vol/RC)*...
                trapz(t, ...
                t.^k.*convPR, 2);
        end

        l = -dkMomIdmua/tkMom;
    end
end