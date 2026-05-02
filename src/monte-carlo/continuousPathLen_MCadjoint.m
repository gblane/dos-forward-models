function [ll, L] = continuousPathLen_MCadjoint(adjoint, t, Vol)
% [ll, L] = continuousPathLen_MCadjoint(adjoint, t, Vol)
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
% 
% Outputs:
%   L       - Total pathlength. (mm)
%   ll      - Matrix of partial pathlengths. (mm)
    
    arguments
        adjoint struct;
        t (:,1) double; %ps
        Vol (1,1) double; %mm^3
    end
    
    PHIsi_CW=trapz(t, adjoint.PHIsi, 4);
    PHIdi_CW=trapz(t, adjoint.PHIdi, 4);
    PHIsd_CW=trapz(t, adjoint.PHIsd, 1);
    
    ll=((PHIsi_CW.*PHIdi_CW)./PHIsd_CW)*Vol;

    L=sum(ll(:));
end