%% Setup
clear; home;

h=6.62607015e-34; %J/Hz
c=299792458; %m/s
AMP=(6.241509074e18)/(1e9); %(charge/s)/nA

lambda=810e-9; %m
Adet=0.565; %mm^2

rho=3; %mm

data(1)=load('fromNEU/parallel-S1-DA_F1.mat');
datNames(1)="par1";
data(2)=load('fromNEU/parallel-S2-DB_F2.mat');
datNames(2)="par2";
data(3)=load('fromNEU/perpendicular-S1-DA_F1.mat');
datNames(3)="per1";
data(4)=load('fromNEU/perpendicular-S2-DB_F2 2.mat');
datNames(4)="per2";

%% Convert to Intensity
Ep=h*c/lambda;
t=data(1).time;
fs=1/median(diff(t));

I=NaN(length(data(1).time), length(data));
InA=NaN(length(data(1).time), length(data));
for i=1:length(data)
    InA(:, i)=-data(i).data;
    I(:, i)=-(data(i).data/data(i).params.igain)*AMP*Ep/Adet; %W/mm^2
end

%% Estimate Noise and Background
Tnoise=1; %sec
Nnoise=round(Tnoise*fs);

noi_all=[];
bck_all=[];
noi_all_nA=[];
bck_all_nA=[];
for i=1:size(I, 2)
    tmp_nA=rmoutliers(InA(:, i));
    tmp=rmoutliers(I(:, i));
    
    bck_all_nA=[bck_all_nA; movmean(tmp_nA, Nnoise)];
    noi_all_nA=[noi_all_nA; movstd(tmp_nA, Nnoise)];

    bck_all=[bck_all; movmean(tmp, Nnoise)];
    noi_all=[noi_all; movstd(tmp, Nnoise)];
end
bck_nA=median(bck_all_nA); %nA
noi_nA=median(noi_all_nA); %nA
noi2bck_nA=noi_nA/bck_nA;
bck=median(bck_all); %W/mm^2
noi=median(noi_all); %W/mm^2
noi2bck=noi/bck;

%% Estimate Peak
pck_all=[];
pck_all_nA=[];
for i=1:size(I, 2)
    tmp=char(datNames(i));
    if strcmp(tmp(1:3), 'per')
        continue;
    end
    
    pck_all_nA=[pck_all_nA; findpeaks(InA(:, i), ...
        'MinPeakHeight', bck_nA+5*noi_nA,...
        'MinPeakDistance', round(1*fs))];
    
    pck_all=[pck_all; findpeaks(I(:, i), ...
        'MinPeakHeight', bck+5*noi,...
        'MinPeakDistance', round(1*fs))];
end
pck=mean(pck_all);
pckMbck=pck-bck;
snr=pckMbck/noi;

pck_nA=mean(pck_all_nA);
pckMbck_nA=pck_nA-bck_nA;
snr_nA=pckMbck_nA/noi_nA;

%% Save
save('Bout_PckNoi.mat');