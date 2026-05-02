function [SD, SS, DS, params] = makeSenMaps_Scat(armt, pert, geo, optProp, opts)
% [SD, SS, DS, params] = makeSenMaps_Scat(armt, pert, geo, optProp, opts)
% 
% Giles Blaney Winter 2020
% 
% NOTE: All units are in mm or 1/mm.
% 
% Coordinate system used by this code assumes z axis (3rd index) is normal
% to the surface of the medium and positive into the medium.
%
% All inputs optional, defaults will be used with no inputs.
% 
% Inputs:   
%   armt    - Arrangement structure with the following feilds:
%               rSrc  - Coordinates of sources (number of sources X 3) {mm}
%                       [x1, y1, z1; x2, y2, z2; ... ; xn, yn, zn]
%               rDet  - Coordinates of detectors (same structure as rSrc)
%                       {mm}
%               SDprs - Single distance pairs (number of pairs X 2)
%                       [sInd1, dInd1; sInd2, dInd2; ... ;sIndn, dIndn]
%               SSprs - Single slope pairs (number of distances X 2 X 
%                       number of pairs)
%                       SSprs(:, :, n)=[sInd1, dInd1; sInd2, dInd2]
%               DSprs - Dual slope pairs (number of distances X 2 X
%                       number of slopes X number of pairs)
%                       SSprs(:, :, 1, n)=[sInd11, dInd11; sInd12, dInd12]
%                       SSprs(:, :, 2, n)=[sInd21, dInd21; sInd22, dInd22]
%   
%   pert    - Perturbation structure with the following feilds:
%               partSz - Size of perturbation volume in voxels (1 X 3)
%                        [nx, ny, nz]
%               dmusp  - Change in the scattering coefficient of 
%                        perturbation {1/mm}
%   
%   geo     - Geometry structure with the following feilds:
%               xMar - Margin in x direction of full volume, volume will
%                      extend from -xMar to +xMar along the x axis {mm}
%               yMar - Margin in y direction of full volume, volume will
%                      extend from -yMar to +yMar along the y axis {mm}
%               zMax - Maximum z of full volume, volume will extend from 0
%                      to +zMax {mm}
%               dr   - Voxel size in three dimentions {mm}
%                      [dx, dy, dz]
%   
%   optProp - Optical properties structure with the following feilds:
%               nin  - Index of refraction inside the medium
%               nout - Index of refraction outside the medium
%               musp - Baseline scattering coefficient within medium {1/mm}
%               mua  - Baseline absorption coefficient within medium {1/mm}
%   
%   opts    - Options structure with the following feilds:
%               fmod - Modulation frequency {Hz}
%               PhN  - Noise in phase measurement {deg}
%               InN  - Noise in intensity measurement {percent}
% 
% Outputs:
%   SD      - Single distnace structure with the following feilds:
%               S_Ph  - Phase sensitivity map (voxels in x X voxels in y X
%                      voxels in z)
%               S_In  - Intensity sensitivity map (same structure as S_Ph)
%               S_PhN - Noise equivalent sensitivity for phase
%               S_InN - Noise equivalent sensitivity for intensity
%   
%   SS      - Single slope structure with the same feilds as SD
%   
%   DS      - Dual slope structure with the same feilds as SD
%   
%   params  - Parameters structure with the following feilds:
%               armt    - Arrangement structure same as input
%               geo     - Geometry structure same as input with the following
%                         additional feilds:
%                           X - Meshgird for x axis (voxels in x X
%                               voxels in y X voxels in z) {mm}
%                           Y - Meshgird for y axis (same structure as X)
%                               {mm}
%                           Z - Meshgird for z axis (same structure as X)
%                               {mm}
%                           x - x axis coordinates (1 X voxels in x) {mm}
%                           y - y axis coordinates (1 X voxels in y) {mm}
%                           z - z axis coordinates (1 X voxels in z) {mm}
%               optProp - Optical properties structure same as input
%               opts    - Options structure same as input

    %% Unpackage Input
    if nargin<=0
        armt.rSrc=[...
            -30, 0, 0;...
            30, 0, 0]; %mm
        armt.rDet=[...
            -5, 0, 0;...
            5, 0, 0]; %mm
        
        %[sInd1, dInd1; ...; sIndn, dIndn]
        armt.SDprs=[...
            1, 2];
        
        %SSprs(:, :, n)=[sInd1, dInd1; sInd2, dInd2]
        armt.SSprs(:, :, 1)=[...
            1, 1;...
            1, 2];
        
        %SSprs(:, :, SSInd, n)=[sInd1, dInd1; sInd2, dInd2]
        armt.DSprs(:, :, 1, 1)=[...
            1, 1;...
            1, 2];
        armt.DSprs(:, :, 2, 1)=[...
            2, 2;...
            2, 1];
    end
    if nargin<=1
        pert.partSz=[20, 20, 4]; %Num Voxels [x, y, z]
        pert.dmusp=0.12; %1/mm
    end
    if nargin<=2
        geo.xMar=45; %mm
        geo.yMar=15; %mm
        geo.zMax=25; %mm
        geo.dr=[0.5, 0.5, 0.5]; %mm
    end
    if nargin<=3
        optProp.nin=1.4;
        optProp.nout=1;
        optProp.musp=1.2; %1/mm
        optProp.mua=0.01; %1/mm

    end
    if nargin<=4
        opts.fmod=1.40625e8; %Hz
        opts.PhN=0.06; %deg
        opts.InN=0.4; %Percent
    end
    
    rSrc_true=armt.rSrc; %mm
    rDet=armt.rDet; %mm
    SDprs=armt.SDprs; %[sInd1, dInd1; ...; sIndn, dIndn]
    SSprs=armt.SSprs; %SSprs(:, :, n)=[sInd1, dInd1; sInd2, dInd2]
    DSprs=armt.DSprs; %SSprs(:, :, SSInd, n)=[sInd1, dInd1; sInd2, dInd2]
    
    xMar=geo.xMar; %mm
    yMar=geo.yMar; %mm
    zMax=geo.zMax; %mm
    dx=geo.dr(1); %mm
    dy=geo.dr(2); %mm
    dz=geo.dr(3); %mm
    
    partSz(1)=pert.partSz(2); %Num Voxels [y, x, z]
    partSz(2)=pert.partSz(1);
    partSz(3)=pert.partSz(3);
    dmusp=pert.dmusp; %1/mm

    mua=optProp.mua; %1/mm
    musp=optProp.musp; %1/mm
    
    omega=2*pi*opts.fmod; %rad/sec
    PhN=opts.PhN*pi/180; %rad
    InN=opts.InN/100; %Frac
    
    rSrc=rSrc_true+[0, 0, 1/musp];

    %% Calculate Pathlengths
    initVar=NaN(size(rSrc, 1), size(rDet, 1));
    L=initVar;
    LIn=initVar;
    LPh=initVar;
    R=initVar;
    for sInd=1:size(rSrc, 1)
        for dInd=1:size(rDet, 1)
            [L(sInd, dInd), R(sInd, dInd)]=...
                complexTotPathLen_Scat(rSrc(sInd, :), rDet(dInd, :), omega, optProp);

            LIn(sInd, dInd)=real(L(sInd, dInd));
            LPh(sInd, dInd)=imag(L(sInd, dInd));
        end
    end

    x=-xMar:dx:xMar;
    y=-yMar:dy:yMar;
    z=0:dz:zMax;

    [X, Y, Z]=meshgrid(x, y, z);
    r=[X(:), Y(:), Z(:)];

    V=dx*dy*dz;

    initVar=NaN(size(r, 1), size(rSrc, 1), size(rDet, 1));
    l=initVar;
    lIn=initVar;
    lPh=initVar;
    for sInd=1:size(rSrc, 1)
        for dInd=1:size(rDet, 1)
            l(:, sInd, dInd)=complexPartPathLen_Scat(...
                rSrc(sInd, :), r, rDet(dInd, :), V, omega, optProp);

            lIn(:, sInd, dInd)=real(l(:, sInd, dInd));
            lPh(:, sInd, dInd)=imag(l(:, sInd, dInd));
        end
    end
    lIn(isnan(lIn))=0;
    lPh(isnan(lPh))=0;

    %% Sen
    % SD
    if ~isempty(SDprs)
        initVar=NaN(size(r, 1), size(SDprs, 1));
        S_SDIn=initVar;
        S_SDPh=initVar;
        for pInd=1:size(SDprs, 1)
            sInd=SDprs(pInd, 1);
            dInd=SDprs(pInd, 2);
            if dInd==0 || sInd==0
                continue;
            end

            S_SDIn(:, pInd)=lIn(:, sInd, dInd)./LIn(sInd, dInd);
            S_SDPh(:, pInd)=lPh(:, sInd, dInd)./LPh(sInd, dInd);
        end
    end

    % SS
    if ~isempty(SSprs)
        initVar=NaN(size(r, 1), size(SSprs, 3));
        S_SSIn=initVar;
        S_SSPh=initVar;
        S_SSIn_num=initVar;
        S_SSPh_num=initVar;
        S_SSIn_den=initVar;
        S_SSPh_den=initVar;
        for pInd=1:size(SSprs, 3)
            sInds=SSprs(:, 1, pInd);
            dInds=SSprs(:, 2, pInd);            
            if sum(dInds==0)>0 || sum(sInds==0)>0
                continue;
            end

            rsd=vecnorm((rDet(dInds, :)-rSrc(sInds, :))');
            rsd_avg=mean(rsd);
            rsd_diff=rsd-rsd_avg;

            S_SSIn_num(:, pInd)=0;
            S_SSPh_num(:, pInd)=0;
            S_SSIn_den(:, pInd)=0;
            S_SSPh_den(:, pInd)=0;
            for rInd=1:length(rsd_diff)        
                S_SSIn_num(:, pInd)=S_SSIn_num(:, pInd)+...
                    (rsd_diff(rInd).*lIn(:, sInds(rInd), dInds(rInd)));
                S_SSIn_den(:, pInd)=S_SSIn_den(:, pInd)+...
                    (rsd_diff(rInd).*LIn(sInds(rInd), dInds(rInd)));

                S_SSPh_num(:, pInd)=S_SSPh_num(:, pInd)+...
                    (rsd_diff(rInd).*lPh(:, sInds(rInd), dInds(rInd)));
                S_SSPh_den(:, pInd)=S_SSPh_den(:, pInd)+...
                    (rsd_diff(rInd).*LPh(sInds(rInd), dInds(rInd)));
            end
            S_SSIn(:, pInd)=S_SSIn_num(:, pInd)./S_SSIn_den(:, pInd);
            S_SSPh(:, pInd)=S_SSPh_num(:, pInd)./S_SSPh_den(:, pInd);
        end
    end

    % DS
    if ~isempty(DSprs)
        initVar=NaN(size(r, 1), size(DSprs, 4));
        S_DSIn=initVar;
        S_DSPh=initVar;
        for SInd=1:size(DSprs, 4)
            DS_SSprs=DSprs(:, :, :, SInd);

            initVar=NaN(size(r, 1), size(DS_SSprs, 3));
            S_DS_SSIn=initVar;
            S_DS_SSPh=initVar;
            S_DS_SSIn_num=initVar;
            S_DS_SSPh_num=initVar;
            S_DS_SSIn_den=initVar;
            S_DS_SSPh_den=initVar;
            for pInd=1:size(DS_SSprs, 3)
                sInds=DS_SSprs(:, 1, pInd);
                dInds=DS_SSprs(:, 2, pInd);
                if sum(dInds==0)>0 || sum(sInds==0)>0
                    continue;
                end

                rsd=vecnorm((rDet(dInds, :)-rSrc(sInds, :))');
                rsd_avg=mean(rsd);
                rsd_diff=rsd-rsd_avg;

                S_DS_SSIn_num(:, pInd)=0;
                S_DS_SSPh_num(:, pInd)=0;
                S_DS_SSIn_den(:, pInd)=0;
                S_DS_SSPh_den(:, pInd)=0;
                for rInd=1:length(rsd_diff)        
                    S_DS_SSIn_num(:, pInd)=S_DS_SSIn_num(:, pInd)+...
                        (rsd_diff(rInd).*lIn(:, sInds(rInd), dInds(rInd)));
                    S_DS_SSIn_den(:, pInd)=S_DS_SSIn_den(:, pInd)+...
                        (rsd_diff(rInd).*LIn(sInds(rInd), dInds(rInd)));

                    S_DS_SSPh_num(:, pInd)=S_DS_SSPh_num(:, pInd)+...
                        (rsd_diff(rInd).*lPh(:, sInds(rInd), dInds(rInd)));
                    S_DS_SSPh_den(:, pInd)=S_DS_SSPh_den(:, pInd)+...
                        (rsd_diff(rInd).*LPh(sInds(rInd), dInds(rInd)));
                end
                S_DS_SSIn(:, pInd)=S_DS_SSIn_num(:, pInd)./S_DS_SSIn_den(:, pInd);
                S_DS_SSPh(:, pInd)=S_DS_SSPh_num(:, pInd)./S_DS_SSPh_den(:, pInd);
            end

            S_DSIn(:, SInd)=nanmean(S_DS_SSIn, 2);
            S_DSPh(:, SInd)=nanmean(S_DS_SSPh, 2);
        end
    end

    %% Make Maps
    % SD
    if ~isempty(SDprs)
        initVar=NaN(size(X, 1), size(X, 2), size(X, 3), size(SDprs, 1));
        S_SDIn_Map=initVar;
        S_SDPh_Map=initVar;
        for pInd=1:size(SDprs, 1)
            S_SDIn_Map(:, :, :, pInd)=reshape(S_SDIn(:, pInd), size(X));
            S_SDPh_Map(:, :, :, pInd)=reshape(S_SDPh(:, pInd), size(X));
        end
    end

    % SS
    if ~isempty(SSprs)
        initVar=NaN(size(X, 1), size(X, 2), size(X, 3), size(SSprs, 3));
        S_SSIn_Map=initVar;
        S_SSPh_Map=initVar;
        for pInd=1:size(SSprs, 3)
            S_SSIn_Map(:, :, :, pInd)=reshape(S_SSIn(:, pInd), size(X));
            S_SSPh_Map(:, :, :, pInd)=reshape(S_SSPh(:, pInd), size(X));
        end
    end

    % DS
    if ~isempty(DSprs)
        initVar=NaN(size(X, 1), size(X, 2), size(X, 3), size(DSprs, 4));
        S_DSIn_Map=initVar;
        S_DSPh_Map=initVar;
        for pInd=1:size(DSprs, 4)
            S_DSIn_Map(:, :, :, pInd)=reshape(S_DSIn(:, pInd), size(X));
            S_DSPh_Map(:, :, :, pInd)=reshape(S_DSPh(:, pInd), size(X));
        end
    end

    %% Convolve Volume
    % SD
    if ~isempty(SDprs)
        initVar=NaN(size(S_SDIn_Map));
        S_SDIn_Map_vol=initVar;
        S_SDPh_Map_vol=initVar;
        for pInd=1:size(SDprs, 1)
            H=ones(partSz);
            S_SDIn_Map_vol(:, :, :, pInd)=convn(...
                S_SDIn_Map(:, :, :, pInd), H, 'same');
            S_SDPh_Map_vol(:, :, :, pInd)=convn(...
                S_SDPh_Map(:, :, :, pInd), H, 'same');
        end
    end

    % SS
    if ~isempty(SSprs)
        initVar=NaN(size(S_SSIn_Map));
        S_SSIn_Map_vol=initVar;
        S_SSPh_Map_vol=initVar;
        for pInd=1:size(SSprs, 3)
            H=ones(partSz);
            S_SSIn_Map_vol(:, :, :, pInd)=convn(...
                S_SSIn_Map(:, :, :, pInd), H, 'same');
            S_SSPh_Map_vol(:, :, :, pInd)=convn(...
                S_SSPh_Map(:, :, :, pInd), H, 'same');
        end
    end

    % DS
    if ~isempty(DSprs)
        initVar=NaN(size(S_DSIn_Map));
        S_DSIn_Map_vol=initVar;
        S_DSPh_Map_vol=initVar;
        for pInd=1:size(DSprs, 4)
            H=ones(partSz);
            S_DSIn_Map_vol(:, :, :, pInd)=convn(...
                S_DSIn_Map(:, :, :, pInd), H, 'same');
            S_DSPh_Map_vol(:, :, :, pInd)=convn(...
                S_DSPh_Map(:, :, :, pInd), H, 'same');
        end
    end
    
    %% Calc Noise
    % SD
    if ~isempty(SDprs)
        initVar=NaN(size(SDprs, 1), 1);
        S_SDInN=initVar;
        S_SDPhN=initVar;
        for pInd=1:size(SDprs, 1)
            sInd=SDprs(pInd, 1);
            dInd=SDprs(pInd, 2);
            if dInd==0 || sInd==0
                continue;
            end

            S_SDInN(pInd)=InN/abs(LIn(sInd, dInd)*dmusp);
            S_SDPhN(pInd)=PhN/abs(LPh(sInd, dInd)*dmusp);
        end
    end
    
    % SS
    if ~isempty(SSprs)
        initVar=NaN(size(SSprs, 3), 1);
        S_SSInN=initVar;
        S_SSPhN=initVar;
        for pInd=1:size(SSprs, 3)
            sInds=SSprs(:, 1, pInd);
            dInds=SSprs(:, 2, pInd);
            if sum(dInds==0)>0 || sum(sInds==0)>0
                continue;
            end
            
            rsd=vecnorm((rDet(dInds, :)-rSrc(sInds, :))');
            rsd_avg=mean(rsd);
            rsd_diff=rsd-rsd_avg;
            
            initVar=NaN(size(rsd_diff));
            LInVec=initVar;
            LPhVec=initVar;
            for rInd=1:length(rsd_diff)
                LInVec(rInd)=LIn(sInds(rInd), dInds(rInd));
                LPhVec(rInd)=LPh(sInds(rInd), dInds(rInd));
            end
            
            S_SSInN(pInd)=(InN*sqrt(sum(rsd_diff.^2)))/...
                (dmusp*dot(rsd_diff, abs(LInVec)));
            S_SSPhN(pInd)=(PhN*sqrt(sum(rsd_diff.^2)))/...
                (dmusp*dot(rsd_diff, abs(LPhVec)));

        end
    end
    
    % DS
    if ~isempty(DSprs)
        initVar=NaN(size(DSprs, 4), 1);
        S_DSInN=initVar;
        S_DSPhN=initVar;
        for SInd=1:size(DSprs, 4)
            DS_SSprs=DSprs(:, :, :, SInd);

            initVar=NaN(size(DS_SSprs, 3), 1);
            S_DS_SSInN=initVar;
            S_DS_SSPhN=initVar;
            for pInd=1:size(DS_SSprs, 3)
                sInds=DS_SSprs(:, 1, pInd);
                dInds=DS_SSprs(:, 2, pInd);
                if sum(dInds==0)>0 || sum(sInds==0)>0
                    continue;
                end

                rsd=vecnorm((rDet(dInds, :)-rSrc(sInds, :))');
                rsd_avg=mean(rsd);
                rsd_diff=rsd-rsd_avg;

                initVar=NaN(size(rsd_diff));
                LInVec=initVar;
                LPhVec=initVar;
                for rInd=1:length(rsd_diff)
                    LInVec(rInd)=LIn(sInds(rInd), dInds(rInd));
                    LPhVec(rInd)=LPh(sInds(rInd), dInds(rInd));
                end

                S_DS_SSInN(pInd)=(InN*sqrt(sum(rsd_diff.^2)))/...
                    (dmusp*dot(rsd_diff, abs(LInVec)));
                S_DS_SSPhN(pInd)=(PhN*sqrt(sum(rsd_diff.^2)))/...
                    (dmusp*dot(rsd_diff, abs(LPhVec)));
            end

            S_DSInN(SInd)=sqrt(nansum(S_DS_SSInN.^2))/length(S_DS_SSInN);
            S_DSPhN(SInd)=sqrt(nansum(S_DS_SSPhN.^2))/length(S_DS_SSPhN);
        end
    end
    
    %% Package Output
    geo.X=X;
    geo.Y=Y;
    geo.Z=Z;
    geo.x=x;
    geo.y=y;
    geo.z=z;
    
    params.armt=armt;
    params.geo=geo;
    params.optProp=optProp;
    params.opts=opts;
    
    if ~isempty(SDprs)
        SD.S_Ph=S_SDPh_Map_vol;
        SD.S_In=S_SDIn_Map_vol;
        SD.S_PhN=S_SDPhN;
        SD.S_InN=S_SDInN;
    else
        SD=[];
    end
    
    if ~isempty(SSprs)
        SS.S_Ph=S_SSPh_Map_vol;
        SS.S_In=S_SSIn_Map_vol;
        SS.S_PhN=S_SSPhN;
        SS.S_InN=S_SSInN;
    else
        SS=[];
    end
    
    if ~isempty(DSprs)
        DS.S_Ph=S_DSPh_Map_vol;
        DS.S_In=S_DSIn_Map_vol;
        DS.S_PhN=S_DSPhN;
        DS.S_InN=S_DSInN;
    else
        DS=[];
    end
end