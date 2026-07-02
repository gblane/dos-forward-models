function [ll, L] = temporalGatePathLen_MCadjoint(adjoint, t, Vol, tg)
% temporalGatePathLen_MCadjoint Calculate temporal gate pathlengths from MC adjoint simulations.
%
% [ll, L] = temporalGatePathLen_MCadjoint(adjoint, t, Vol, tg)
%
% Written by Giles Blaney, Ph.D. Summer 2023
%
% Inputs:
%   adjoint - Adjoint simulation results structure [struct]
%   t       - Time vector [ps]
%   Vol     - Voxel volume [mm^3]
%   tg      - Gate start and end time [ps]
%
% Outputs:
%   ll - Matrix of partial pathlengths [mm]
%   L  - Total pathlength [mm]

    arguments
        adjoint struct;
        t (:,1) double; % ps
        Vol (1,1) double; % mm^3
        tg (1,2) double; % ps
    end

    dt = median(diff(t));

    convPR = NaN(size(adjoint.PHIsi, 1), size(adjoint.PHIsi, 2), ...
        size(adjoint.PHIsi, 3), length(t));
    for i = 1:size(adjoint.PHIsi, 1)
        for j = 1:size(adjoint.PHIsi, 2)
            for k = 1:size(adjoint.PHIsi, 3)
                convPR(i, j, k, :) = ...
                    conv(...
                    squeeze(adjoint.PHIsi(i, j, k, :)),...
                    squeeze(adjoint.PHIdi(i, j, k, :)),...
                    "same")*dt;
            end
        end
    end

    [~, i1] = min(abs(t-tg(1)));
    [~, i2] = min(abs(t-tg(2)));

    PHIsd_g = trapz(t(i1:i2), adjoint.PHIsd(i1:i2));
    num = trapz(t(i1:i2), convPR(:, :, :, i1:i2), 4);
    ll = (num./PHIsd_g)*Vol;

    L = sum(ll(:));
end
