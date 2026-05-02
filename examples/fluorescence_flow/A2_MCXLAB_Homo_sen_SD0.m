%% Setup
clear; home;

np=1e9;
% np=1e8;

grdSp=0.1; %mm

% Dims
Lx=30+grdSp; %mm (x side length)
Ly=30; %mm (y side length)
Lz=30; %mm (z side length)

r_src=[0, 0, 0]; %mm [x, y, z]
r_det=[0, 0, 0]; %mm [x, y, z]

% Opt Prop for 800 nm
mua=0.002; %1/mm
musp=0.7; %1/mm
g=0.9;
mus=musp/(1-g); %1/mm
nin=1.37;
nout=1;

% Time
ndt=10;
tend=10e-9;

%% MCX Setup
% Req
cfg0.nphoton=np;
cfg0.seed=round(1000*rand);
cfg0.vol=uint8(ones(Lx/grdSp, Ly/grdSp, Lz/grdSp)); %Media ind (0=outside)
cfg0.prop=[...
    0 0 1 nout;... Outside
    mua mus g nin... Inside
    ]; %mua mus g n

cfg0.tstart=0; %sec
cfg0.tstep=tend/ndt;
cfg0.tend=tend;

%srcpos and srcdir defined later

%Opt MC Sim
cfg0.isreflect=nin~=nout; %Consider n mismatch (0 = matched)
cfg0.unitinmm=grdSp; %mm

%Opt GPU
cfg0.autopilot=1;
cfg0.gpuid=1;
cfg0.isgpuinfo=1;

%Opt Output
cfg0.outputtype='flux';
cfg0.issaveref=0;

%Opt Debug
cfg0.debuglevel='p';

%% Make Coor Sys
gridOrigin=[Lx/2, Ly/2, -grdSp/2]+grdSp/2; %grid (mm)

[YY, XX, ZZ]=meshgrid(...
    (grdSp:grdSp:(size(cfg0.vol, 2)*grdSp))-gridOrigin(2),...
    (grdSp:grdSp:(size(cfg0.vol, 1)*grdSp))-gridOrigin(1),...
    (0:grdSp:((size(cfg0.vol, 3)-1)*grdSp))-gridOrigin(3));

x=squeeze(XX(:, 1, 1));
y=squeeze(YY(1, :, 1));
z=squeeze(ZZ(1, 1, :));

t=((cfg0.tstart+cfg0.tstep/2):cfg0.tstep:(cfg0.tend-cfg0.tstep/2))';

%% Place Optodes
cfg0.srcpos=[];
cfg0.srcdir=[];
cfg0.srctype=[];
cfg0.srcparam1=0;

%Det
% cfg(1)=cfg0;
% cfg(1).srcpos=(r_det(1, :)+gridOrigin)/grdSp+[0.5, 0.5, 1];
% cfg(1).srcdir=[0, 0, 1]; %vec
% cfg(1).srctype='isotropic';

%Det
cfg(1)=cfg0;
cfg(1).srcpos=(r_det(1, :)+gridOrigin)/grdSp+[0.5, 0.5, 1];
cfg(1).srcdir=[0, 0, 1]; %vec
cfg(1).srctype='cone';
cfg(1).srcparam1=asin(0.5);

%Src
cfg(2)=cfg0;
cfg(2).srcpos=(r_src(1, :)+gridOrigin)/grdSp+[0.5, 0.5, 1];
cfg(2).srcdir=[0, 0, 1]; %vec
cfg(2).srctype='pencil';

%% Run MCX
tic;
MCout=mcxlab(cfg);
runTime=toc;

%% Calc PHI
% Det, Src
PHI=NaN(length(x), length(y), length(z), 2);
for i=1:length(cfg)
    PHI(:, :, :, i)=trapz(t, MCout(i).data, 4);
end

%% Save
save('MCout_SD0.mat', '-v7.3');