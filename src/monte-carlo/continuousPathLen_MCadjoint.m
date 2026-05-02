function [ll, L] = continuousPathLen_MCadjoint(adjoint, t, Vol)
% continuousPathLen_MCadjoint Calculate CW pathlengths from MC adjoint simulations.
%
% [ll, L] = continuousPathLen_MCadjoint(adjoint, t, Vol)
%
% Written by Giles Blaney, Ph.D. Summer 2023
%
% Inputs:
%   adjoint - Adjoint simulation results structure [struct]
%   t       - Time vector [ps]
%   Vol     - Voxel volume [mm^3]
%
% Outputs:
%   ll - Matrix of partial pathlengths [mm]
%   L  - Total pathlength [mm]
    
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