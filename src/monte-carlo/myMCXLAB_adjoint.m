function [adjoint, params, MCXout] = myMCXLAB_adjoint(rs, rd, optProp, NVA)
% [adjoint, params, MCXout] = myMCXLAB_adjoint(rs, rd, optProp, NVA)
% 
% Giles Blaney Ph.D. Summer 2023
% 
% Inputs:
%   rs       - Pencil beam source coordinates [x, y, z]. (mm)
%   rd       - Detector coordinates in [x, y, z]. (mm)
%   optProp  - (OPTIONAL) Struct of optical properties with the following
%                fields:
%                    nin  - (default=1.333) Index of refraction inside. (-)
%                    nout - (default=1) Index of refraction outside. (-)
%                    musp - (default=1.1 1/mm) Reduced scattering. (1/mm)
%                    g    - (default=0.9) anisotropy. (-)
%                    mua  - (default=0.011 1/mm) Absorption. (1/mm)
% 
% Name Value Arguments:
%   'np'     - (default=1e8) Number of photons.
%   'dr'     - (default=1 mm) Voxel side length. (mm)
%   'xl'     - (default=[-15.5, 15.5] mm) Limits in the x direction in 
%               [min, max]. (mm)
%   'yl'     - (default=[-15.5, 15.5] mm) Limits in the y direction in
%               [min, max]. (mm)
%   'zl'     - (default=[0, 30] mm) Limits in the y direction in 
%               [min, max]. (mm)
%   'ndt'    - (default=1e3) Number of time bins.
%   'tend'   - (default=5e3 ps) Max time bin edge. (ps) 
%   'detNA'  - (default=0.5) Detector NA.
% 
% Outputs:
%   adjoint  - Struct with the following fields:
%               PHIsd   - The fluence from the source to the voxel at the
%                           detector. (1/(ps mm^2))
%               PHIsi   - The fluence from the source to each voxel.
%                           (1/(ps mm^2))
%               PHIdi   - The fluence from the detector (for adjoint) to 
%                           each voxel. (1/(ps mm^2))
%   params   - Struct with the following fields:
%               t       - Time vector. (ps)
%               x       - X-axis vector. (mm)
%               y       - Y-axis vector. (mm)
%               z       - Z-axis vector. (mm)
%               runTime - Time to run MCX. (s)
%   MCXout   - Output struct from mcxlab()
    
    arguments
        rs (1,3) double; %mm [x, y, z]
        rd (1,3) double; %mm [x, y, z]

        optProp struct = [];

        NVA.np (1,1) double = 1e8;

        NVA.dr (1,1) double = 1; %mm
        NVA.xl (1,2) double = [-15.5, 15.5]; %mm
        NVA.yl (1,2) double = [-15.5, 15.5]; %mm
        NVA.zl (1,2) double = [0, 30]; %mm

        NVA.ndt  (1,1) double = 1e3; 
        NVA.tend (1,1) double = 5e3; %ps

        NVA.detNA (1,1) double = 0.5; 
    end

    %% Check Inputs
    % Default optical properties if needed
    if isempty(optProp)
        clear optProp;

        optProp.nin=1.333;
        optProp.nout=1;
        optProp.musp=1.1; %1/mm
        optProp.g=0.9;
        optProp.mua=0.011; %1/mm
        
        warning('Default optical properties used');
    end
    
    if ~exist('mcxlab', 'file')
        error(['mcxlab not found,' ...
            ' MCX must be installed and in path; ' ...
            'see: http://mcx.space/wiki/index.cgi?Doc/MCXLAB']);
    end

    %% Setup MCX cfg0
    % Req
    cfg0.nphoton=NVA.np;
    cfg0.seed=round(1000*rand);
    cfg0.vol=uint8(ones( ...
        (NVA.xl(2)-NVA.xl(1))/NVA.dr, ...
        (NVA.yl(2)-NVA.yl(1))/NVA.dr, ...
        (NVA.zl(2)-NVA.zl(1))/NVA.dr)); %Media ind (0=outside)
    cfg0.prop=[...
                  0,                          0,         1, optProp.nout;...
        optProp.mua, optProp.musp/(1-optProp.g), optProp.g,  optProp.nin ...
        ]; %mua mus g n
    
    cfg0.tstart=0; %s
    cfg0.tstep=(NVA.tend/NVA.ndt)*1e-12; %s
    cfg0.tend=NVA.tend*1e-12; %s
    
    cfg0.srcpos=[];
    cfg0.srcdir=[];
    cfg0.srctype=[];
    cfg0.srcparam1=0;
    
    %Opt MC Sim
    cfg0.isreflect=(optProp.nin~=optProp.nout); %Consider n mismatch
    cfg0.unitinmm=NVA.dr; %mm
    
    %Opt GPU
    cfg0.autopilot=1;
    cfg0.gpuid=1;
    cfg0.isgpuinfo=1;
    
    %Opt Output
    cfg0.outputtype='flux';
    cfg0.issaveref=0;
    
    %Opt Debug
    cfg0.debuglevel='p';
    
    %% Make coordinate systems
    tNeg_sec=(-(cfg0.tend-cfg0.tstep/2): ...
        cfg0.tstep: ...
        (cfg0.tstart-cfg0.tstep/2))';
    tPos_sec=((cfg0.tstart+cfg0.tstep/2): ...
        cfg0.tstep: ...
        (cfg0.tend-cfg0.tstep/2))';
    params.t=[tNeg_sec; tPos_sec]*1e12; %ps
    
    gridOrigin=[ ...
        -NVA.xl(1), ...
        -NVA.yl(1), ...
        -NVA.dr/2]+NVA.dr/2; %grid mm

    [params.YY, params.XX, params.ZZ]=meshgrid(...
        (NVA.dr:NVA.dr:(size(cfg0.vol, 2)*NVA.dr))-gridOrigin(2),...
        (NVA.dr:NVA.dr:(size(cfg0.vol, 1)*NVA.dr))-gridOrigin(1),...
        (0:NVA.dr:((size(cfg0.vol, 3)-1)*NVA.dr))-gridOrigin(3));

    params.x=squeeze(params.XX(:, 1, 1));
    params.y=squeeze(params.YY(1, :, 1));
    params.z=squeeze(params.ZZ(1, 1, :));

    %% Add sources and detectors and make cfg
    %Src
    cfg(1)=cfg0;
    cfg(1).srcpos=(rs+gridOrigin)/NVA.dr+[0.5, 0.5, 1];
    cfg(1).srcdir=[0, 0, 1]; %vec
    cfg(1).srctype='pencil';

    %Det
    cfg(2)=cfg0;
    cfg(2).srcpos=(rd+gridOrigin)/NVA.dr+[0.5, 0.5, 1];
    cfg(2).srcdir=[0, 0, 1]; %vec
    cfg(2).srctype='cone';
    cfg(2).srcparam1=asin(NVA.detNA/optProp.nin);

    %% Run MCX
    tic;
    MCXout=mcxlab(cfg);
    params.runTime=toc;

    %% Parse MCXout
    % MC units in 1/(s mm^2)
    % Convert to 1/(ps mm^2)
    adjoint.PHIsi=zeros(length(params.x), length(params.y), ...
        length(params.z), length(params.t));
    adjoint.PHIdi=zeros(length(params.x), length(params.y), ...
        length(params.z), length(params.t));
    adjoint.PHIsi(:, :, :, params.t>=tPos_sec(1))=MCXout(1).data*1e-12;
    adjoint.PHIdi(:, :, :, params.t>=tPos_sec(1))=MCXout(2).data*1e-12;

    adjoint.PHIsd=squeeze( ...
        adjoint.PHIsi( ...
        round(cfg(2).srcpos(1)-0.5), ...
        round(cfg(2).srcpos(2)-0.5), ...
        cfg(2).srcpos(3), :));

end