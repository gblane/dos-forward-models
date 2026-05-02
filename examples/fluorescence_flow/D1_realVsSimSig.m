%% Setup
clear; home;

SD1Apar=load('fromNEU/parallel-S1-DA_F1.mat');
SD2Bpar=load('fromNEU/parallel-S2-DB_F2.mat');
SD1Aper=load('fromNEU/perpendicular-S1-DA_F1.mat');
SD2Bper=load('fromNEU/perpendicular-S2-DB_F2 2.mat');

load('MCout.mat', 'x', 'y', 'z', 'PHI', 'grdSp');
load('Bout_PckNoi.mat');
load('Bout_hom_eta.mat');

dx=1; %sec
yl=[-0.2, 1.5];
xl=[-dx/2, dx/2]; %sec
vel=25; %mm/sec
Nsm=120;

V=grdSp^3;

z_slice=1.5;

Psrc=75e-3; %W

%% Calc W
W1A=PHI(:, :, :, 2).*PHI(:, :, :, 1)*V;

%% R [W/mm^2]
[~, xInd]=min(abs(x-0));
[~, yInd]=min(abs(y-0));
[~, zInd]=min(abs(z-z_slice));

R1Apck_per=Psrc*W1A(xInd, :, zInd)*etaMuaPck;
R1Apck_par=Psrc*W1A(:, yInd, zInd)*etaMuaPck;

%% Plot
t=SD1Apar.time-75.423;
inds=or(t<=xl(1), t>=xl(2));
t(inds)=[];
I1A=-SD1Apar.data;
I1A(inds)=[];
I1A0=mean(I1A([1:100, (end-100):end]));
tmp=(I1A-I1A0);
SigLP=movmean(tmp, Nsm);
SigLP=SigLP-mean(SigLP([1:100, (end-100):end]));
a1=max(SigLP);

t=SD2Bpar.time-658.679;
inds=or(t<=xl(1), t>=xl(2));
t(inds)=[];
I1A=-SD2Bpar.data;
I1A(inds)=[];
I1A0=mean(I1A([1:100, (end-100):end]));
tmp=(I1A-I1A0);
SigLP=movmean(tmp, Nsm);
SigLP=SigLP-mean(SigLP([1:100, (end-100):end]));
a2=max(SigLP);

a=mean([a1, a2]);

dx=1; %sec
yl=[-0.2, 1.5];
xl=[-dx/2, dx/2]; %sec
vel=25; %mm/sec
Nsm=120;

h=figure(1000); clf;
h.Name='sigVsim';

subaxis(2, 2, 3);
t=SD1Apar.time-75.423;
inds=or(t<=xl(1), t>=xl(2));
t(inds)=[];
I1A=-SD1Apar.data;
I1A(inds)=[];
I1A0=mean(I1A([1:100, (end-100):end]));
tmp=(I1A-I1A0);
SigLP=movmean(tmp, Nsm);
SigLP=SigLP-mean(SigLP([1:100, (end-100):end]));
ax(1)=plot_avg(t, SigLP/a, ...
    ones(size(t))*std(tmp([1:100, (end-100):end])/a),...
    'b', 'b'); hold on;
tmp=movmean(R1Apck_par, ...
    round(Nsm*median(diff(t))/(median(diff(x/vel)))));
b=max(tmp);
ax(2)=plot(x/vel, tmp/b, '--r'); hold off;
xlim(xl);
ylim(yl);
xlabel('$t$ (sec)');
ylabel(sprintf('Normalized\nSingle-Distance (SD)'));
title(sprintf('\\textbf{(e)} Source 1 - Detector A\nParallel Flow'));

subaxis(2, 2, 4);
t=SD2Bpar.time-658.679;
inds=or(t<=xl(1), t>=xl(2));
t(inds)=[];
I1A=-SD2Bpar.data;
I1A(inds)=[];
I1A0=mean(I1A([1:100, (end-100):end]));
tmp=(I1A-I1A0);
SigLP=movmean(tmp, Nsm);
SigLP=SigLP-mean(SigLP([1:100, (end-100):end]));
ax(1)=plot_avg(t, SigLP/a2, ...
    ones(size(t))*std(tmp([1:100, (end-100):end])/a),...
    'b', 'b'); hold on;
tmp=movmean(flipud(R1Apck_par), ...
    round(Nsm*median(diff(t))/(median(diff(x/vel)))));
ax(2)=plot(x/vel, tmp/b, '--r'); hold off;
xlim(xl);
ylim(yl);
xlabel('$t$ (sec)');
set(gca, 'YTickLabel', {});
leg=legend(ax, 'Experimental Data', 'Monte-Carlo Model',...
    'Location', 'northwest', 'Orientation','horizontal');
leg.Position=leg.Position+[-0.21, 0.02, 0, 0];
title(sprintf('\\textbf{(f)} Source 2 - Detector B\nParallel Flow'));

subaxis(2, 2, 1, 'sh', 0.05, 'sv', 0.1,...
    'mb', 0.1, 'mt', 0.1', 'ml', 0.1, 'mr', 0.03);
t=SD1Aper.time-276.44;
inds=or(t<=xl(1), t>=xl(2));
t(inds)=[];
I1A=-SD1Aper.data;
I1A(inds)=[];
I1A0=mean(I1A([1:100, (end-100):end]));
tmp=(I1A-I1A0);
SigLP=movmean(tmp, Nsm);
SigLP=SigLP-mean(SigLP([1:100, (end-100):end]));
ax(1)=plot_avg(t, SigLP/a, ...
    ones(size(t))*std(tmp([1:100, (end-100):end])/a),...
    'b', 'b'); hold on;
tmp=movmean(R1Apck_per, ...
    round(Nsm*median(diff(t))/(median(diff(x/vel)))));
ax(2)=plot(y/vel, tmp/b, '--r'); hold off;
xlim(xl);
ylim(yl);
set(gca, 'XTickLabel', {});
ylabel(sprintf('Normalized\nSingle-Distance (SD)'));
% legend(ax, 'Experimental Data', 'Monte-Carlo Model');
title(sprintf('\\textbf{(b)} Source 1 - Detector A\nPerpendicular Flow'));

subaxis(2, 2, 2);
t=SD2Bper.time-425.465;
inds=or(t<=xl(1), t>=xl(2));
t(inds)=[];
I1A=-SD2Bper.data;
I1A(inds)=[];
I1A0=mean(I1A([1:100, (end-100):end]));
tmp=(I1A-I1A0);
SigLP=movmean(tmp, Nsm);
SigLP=SigLP-mean(SigLP([1:100, (end-100):end]));
ax(1)=plot_avg(t, SigLP/a, ...
    ones(size(t))*std(tmp([1:100, (end-100):end])/a),...
    'b', 'b'); hold on;
tmp=movmean(R1Apck_per, ...
    round(Nsm*median(diff(t))/(median(diff(x/vel)))));
ax(2)=plot(y/vel, tmp/b, '--r'); hold off;
xlim(xl);
ylim(yl);
set(gca, 'XTickLabel', {});
set(gca, 'YTickLabel', {});
title(sprintf('\\textbf{(c)} Source 2 - Detector B\nPerpendicular Flow'));