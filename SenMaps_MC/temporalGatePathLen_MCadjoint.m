function [ll, L] = temporalGatePathLen_MCadjoint(adjoint, t, Vol, tg)
% [ll, L] = temporalGatePathLen_MCadjoint(adjoint, t, Vol, tg)
% 
% Giles Blaney Ph.D. Summer 2023
% 
% Inputs:
%   adjoint  - Struct with the following fields:
%               PHIsd   - The fluence from the source to the voxel at the
%                           detector. (1/(ps mm^2))
%               PHIsi   - The fluence from the source to each voxel.
%                           (1/(ps mm^2))
%               PHIdi   - The fluence from the detector (for adjoint) to 
%                           each voxel. (1/(ps mm^2))
%   t        - Time. (ps)
%   Vol      - Volume. (mm^3)
%   tg       - Gate start and end time. (ps)
% 
% Outputs:
%   L        - Total pathlength. (mm)
%   ll       - Matrix of partial pathlengths. (mm)
    
    arguments
        adjoint struct;
        t (:,1) double; %ps
        Vol (1,1) double; %mm^3
        tg (1,2) double; %ps
    end

    dt=median(diff(t));
    
    convPR=NaN(size(adjoint.PHIsi, 1), size(adjoint.PHIsi, 2), ...
        size(adjoint.PHIsi, 3), length(t));
    for i=1:size(adjoint.PHIsi, 1)
        for j=1:size(adjoint.PHIsi, 2)
            for k=1:size(adjoint.PHIsi, 3)
                convPR(i, j, k, :)=...
                    conv(...
                    squeeze(adjoint.PHIsi(i, j, k, :)),...
                    squeeze(adjoint.PHIdi(i, j, k, :)),...
                    'same')*dt;
            end
        end
    end

    [~, i1]=min(abs(t-tg(1)));
    [~, i2]=min(abs(t-tg(2)));
    
    PHIsd_g=trapz(t(i1:i2), adjoint.PHIsd(i1:i2));
    num=trapz(t(i1:i2), convPR(:, :, :, i1:i2), 4);
    ll=(num./PHIsd_g)*Vol;

    L=sum(ll(:));
end