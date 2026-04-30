function [SD, SS, DS, params] = simMeasFromSenMap(simParams, senMapParams)
% Giles Blaney Fall 2019

    %% Unpackage Input
    if nargin<=0
        %[dmua_1_1, ..., dmua_n_1; ...; dmua_1_tn, ..., dmua_n_tn]
        dmua=[0, 0;...
            0.0015, -0.0015;...
            0.003, -0.003]; %1/mm
        
        %bnds(:, :, n, tn)=[x_min, y_min, z_min; x_max, y_max, z_max]
        bndsTemp(:, :, 1)=[-Inf, -Inf, 0; Inf, Inf, 10]; %mm
        bndsTemp(:, :, 2)=[-Inf, -Inf, 15; Inf, Inf, Inf]; %mm
        bnds=repmat(bndsTemp, 1, 1, 1, size(dmua, 1));
    else
        bnds=simParams.bnds;
        dmua=simParams.dmua;
    end
    if nargin<=1
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

        pert.partSz=[1, 1, 1]; %Num Voxels [x, y, z]
        pert.dmua=0.003*(20*20*4); %1/mm

        geo.xMar=45; %mm
        geo.yMar=45; %mm
        geo.zMax=25; %mm
        geo.dr=[0.5, 0.5, 0.5]; %mm

        optProp.nin=1.4;
        optProp.nout=1;
        optProp.musp=1.2; %1/mm
        optProp.mua=0.01; %1/mm
        
        opts.fmod=1.40625e8; %Hz
        opts.PhN=0.06; %deg
        opts.InN=0.4; %Percent
    else
        armt=senMapParams.armt;
        pert=senMapParams.pert;
        geo=senMapParams.geo;
        optProp=senMapParams.optProp;
        opts=senMapParams.opts;
    end
    
    pert.partSz=[1, 1, 1];
    
    %% Make Sen Maps
    [SD, SS, DS, params]=makeSenMaps(armt, pert, geo, optProp, opts);
    X=params.geo.X;
    Y=params.geo.Y;
    Z=params.geo.Z;
    
    %% Make dmua Map
    dmuaMaps=zeros([size(X), size(dmua, 1)]);
    for tn=1:size(dmuaMaps, 4)
        for n=1:size(dmua, 2)
            Xinds=and(X>=bnds(1, 1, n, tn), X<=bnds(2, 1, n, tn));
            Yinds=and(Y>=bnds(1, 2, n, tn), Y<=bnds(2, 2, n, tn));
            Zinds=and(Z>=bnds(1, 3, n, tn), Z<=bnds(2, 3, n, tn));
            inds=and(and(Xinds, Yinds), Zinds);
            
            dmuaMapTemp=zeros(size(X));
            dmuaMapTemp(inds)=dmua(tn, n);
            
            dmuaMaps(:, :, :, tn)=dmuaMaps(:, :, :, tn)+dmuaMapTemp;
        end
    end
    
    %% Find dmua Measured
    if ~isempty(SD)
        for tn=1:size(dmuaMaps, 4)
            dmuaMeasTemp=SD.S_In.*dmuaMaps(:, :, :, tn);
            SD.dmuaMeas_In(tn)=sum(dmuaMeasTemp(:));
            
            dmuaMeasTemp=SD.S_Ph.*dmuaMaps(:, :, :, tn);
            SD.dmuaMeas_Ph(tn)=sum(dmuaMeasTemp(:));
        end
    end
    
    if ~isempty(SS)
        for tn=1:size(dmuaMaps, 4)
            dmuaMeasTemp=SS.S_In.*dmuaMaps(:, :, :, tn);
            SS.dmuaMeas_In(tn)=sum(dmuaMeasTemp(:));
            
            dmuaMeasTemp=SS.S_Ph.*dmuaMaps(:, :, :, tn);
            SS.dmuaMeas_Ph(tn)=sum(dmuaMeasTemp(:));
        end
    end
    
    if ~isempty(DS)
        for tn=1:size(dmuaMaps, 4)
            dmuaMeasTemp=DS.S_In.*dmuaMaps(:, :, :, tn);
            DS.dmuaMeas_In(tn)=sum(dmuaMeasTemp(:));
            
            dmuaMeasTemp=DS.S_Ph.*dmuaMaps(:, :, :, tn);
            DS.dmuaMeas_Ph(tn)=sum(dmuaMeasTemp(:));
        end
    end
    
    %% Package Output
    params.sim.bnds=bnds;
    params.sim.dmua=dmua;
    params.sim.dmuaMaps=dmuaMaps;
end

