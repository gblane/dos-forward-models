function [SD, SS, DS, params] = makeSenMaps2L(armt, pert, geo, optProp, opts)
% makeSenMaps2L Generate sensitivity maps for a two-layer medium.
%
% [SD, SS, DS, params] = makeSenMaps2L(armt, pert, geo, optProp, opts)
%
% Written by Giles Blaney (Spring 2020; Ph.D. awarded May 2022)
%
% Inputs:
%   armt    - Arrangement structure [struct]
%   pert    - Perturbation structure [struct]
%   geo     - Geometry structure [struct]
%   optProp - Optical properties structure [struct]
%   opts    - Options structure [struct]
%
% Outputs:
%   SD     - Single distance sensitivity structure [struct]
%   SS     - Single slope sensitivity structure [struct]
%   DS     - Dual slope sensitivity structure [struct]
%   params - Parameters structure [struct]

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
        pert.partSz=[1, 1, 1]; %Num Voxels [x, y, z]
%         pert.partSz=[20, 20, 1]; %Num Voxels [x, y, z]
        pert.dmua=0.010; %1/mm
    end
    if nargin<=2
        geo.xMar=45; %mm
        geo.yMar=30; %mm
        geo.zMin=0.5; %mm
        geo.zMax=19.5; %mm
        geo.thk=10; %mm
        geo.dr=[1, 1, 1]; %mm
    end
    if nargin<=3
        optProp.nin=[1.4, 1.4];
        optProp.nout=1;
        optProp.musp=[1.20, 0.25]; %1/mm
        optProp.mua=[0.008, 0.020]; %1/mm
    end
    if nargin<=4
        opts.fmod=1.40625e8; %Hz
        opts.PhN=0.04; %deg
        opts.InN=0.4; %Percent
        opts.h_end=500;
        opts.B=50; %mm
        opts.totL='2L';
    end
    
    en=zeroOrdBesselRoots(opts.h_end);
    
    rSrc_true=armt.rSrc; %mm
    rDet=armt.rDet; %mm
    SDprs=armt.SDprs; %[sInd1, dInd1; ...; sIndn, dIndn]
    SSprs=armt.SSprs; %SSprs(:, :, n)=[sInd1, dInd1; sInd2, dInd2]
    DSprs=armt.DSprs; %SSprs(:, :, SSInd, n)=[sInd1, dInd1; sInd2, dInd2]
    
    xMar=geo.xMar; %mm
    yMar=geo.yMar; %mm
    zMin=geo.zMin; %mm
    zMax=geo.zMax; %mm
    dx=geo.dr(1); %mm
    dy=geo.dr(2); %mm
    dz=geo.dr(3); %mm
    
    partSz(1)=pert.partSz(2); %Num Voxels [y, x, z]
    partSz(2)=pert.partSz(1);
    partSz(3)=pert.partSz(3);
    dmua=pert.dmua; %1/mm

    musp=optProp.musp; %1/mm
    
    PhN=opts.PhN*pi/180; %rad
    InN=opts.InN/100; %Frac
    
    rSrc=rSrc_true+[0, 0, 1/musp(1)];

    %% Calculate Pathlengths
    if strcmp(opts.totL, 'homoEff')
        i=0;
        rhosTmp_all=NaN(1, size(rSrc, 1)*size(rDet, 1));
        for sInd=1:size(rSrc, 1)
            for dInd=1:size(rDet, 1)
                i=i+1;
                rhosTmp_all(i)=vecnorm(rDet(dInd, :)-rSrc(sInd, :));
            end
        end
        rhosTmp=unique(rhosTmp_all);
        
        warning('Assuming homogenous nin for total L');
        
        optsEffHomo.fmod=opts.fmod; %Hz
        optsEffHomo.ni=mean(optProp.nin);
        optsEffHomo.no=optProp.nout;
        optsEffHomo.B=200; %mm
        optsEffHomo.h_end=3000;
        homoY=twoLayEffHomoOptProp([optProp.mua, optProp.musp, geo.thk],...
            rhosTmp, optsEffHomo);

        optPropHomo.nin=optsEffHomo.ni;
        optPropHomo.nout=optsEffHomo.no;
        optPropHomo.mua=homoY(1);
        optPropHomo.musp=homoY(2);

        optProp.muaHomo=optPropHomo.mua;
        optProp.muspHomo=optPropHomo.musp;
    end
    
    initVar=NaN(size(rSrc, 1), size(rDet, 1));
    L=initVar;
    LIn=initVar;
    LPh=initVar;
    for sInd=1:size(rSrc, 1)
        for dInd=1:size(rDet, 1)
            switch opts.totL
                case '2L'
                    [L(sInd, dInd), ~, ~]=...
                        complexTotPathLen2L(rSrc(sInd, :), rDet(dInd, :),...
                        geo.thk, en, optProp, opts);
                case 'homoEff'
                    [L(sInd, dInd), ~]=...
                        complexTotPathLen(rSrc(sInd, :), rDet(dInd, :),...
                        2*pi*opts.fmod, optPropHomo);
                otherwise
                    error('Unknown opts.totL');
            end
            
            LIn(sInd, dInd)=real(L(sInd, dInd));
            LPh(sInd, dInd)=imag(L(sInd, dInd));
        end
    end

    x=-xMar:dx:xMar;
    y=-yMar:dy:yMar;
    z=zMin:dz:zMax;

    [X, Y, Z]=meshgrid(x, y, z);
    r=[X(:), Y(:), Z(:)];

    V=dx*dy*dz;

    initVar=NaN(size(r, 1), size(rSrc, 1), size(rDet, 1));
    l=initVar;
    lIn=initVar;
    lPh=initVar;
    for sInd=1:size(rSrc, 1)
        for dInd=1:size(rDet, 1)
            l(:, sInd, dInd)=complexPartPathLen2L(...
                rSrc(sInd, :), r, rDet(dInd, :),...
                V, geo.thk, en, optProp, opts);
            
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

            S_SDInN(pInd)=InN/(LIn(sInd, dInd)*dmua);
            S_SDPhN(pInd)=PhN/(LPh(sInd, dInd)*dmua);
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
                (dmua*dot(rsd_diff, LInVec));
            S_SSPhN(pInd)=(PhN*sqrt(sum(rsd_diff.^2)))/...
                (dmua*dot(rsd_diff, LPhVec));

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
                    (dmua*dot(rsd_diff, LInVec));
                S_DS_SSPhN(pInd)=(PhN*sqrt(sum(rsd_diff.^2)))/...
                    (dmua*dot(rsd_diff, LPhVec));
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
