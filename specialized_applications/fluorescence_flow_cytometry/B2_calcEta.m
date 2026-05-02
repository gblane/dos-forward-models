%% Setup
home; clear;

load('MCout.mat', 'x', 'y', 'z', 'PHI', 'ZZ');
load('Bout_PckNoi.mat');

V=median(diff(x))*median(diff(y))*median(diff(z)); %mm^3

z_bead=1.5; %mm
Psrc=75e-3; %W

%% Make W
% rho=3 mm
W3=PHI(:, :, :, 2).*PHI(:, :, :, 1)*V; %1/mm

%% Find Geometric Inds and Peak W
[~, y0Ind]=min(abs(y-0));
[~, zBeadInd]=min(abs(z-z_bead));

% Max W loc
[~, xPkInd]=max(W3(:, y0Ind, zBeadInd));

%% Sim Loop
for SimNm=["hom", "het"]
    switch SimNm
        case "hom"
            afDist=ones(size(ZZ));
        case "het"
            afDist=exp(log(0.5)*ZZ/0.1);
        otherwise
            afDist=[];
    end
    afDist=afDist/max(afDist(:));
    
    %% Calc Eta
    % Background
    etaMuaBck=(bck/Psrc)/(sum(afDist(:).*W3(:))); %1/mm
    
    % Peak
    etaMuaPck=(pckMbck/Psrc)/W3(xPkInd, y0Ind, zBeadInd); %1/mm
    
    %% Save
    save(join(['Bout_' SimNm '_eta.mat'], ''), ...
        'etaMuaBck', 'etaMuaPck', 'afDist');
end