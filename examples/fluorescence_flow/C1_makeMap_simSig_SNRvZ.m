%% Setup
home; clear;

load('MCout_SD0.mat', 'x', 'y', 'z', 'PHI', 'grdSp');
PHI_SD0=PHI;
x_SD0=x; y_SD0=y; z_SD0=z;
grdSp_SD0=grdSp;
load('MCout.mat', 'x', 'y', 'z', 'PHI', 'grdSp');

load('Bout_PckNoi.mat');
load('Bout_het_eta.mat');
etaMuaBck_het=etaMuaBck;
afDist_het=afDist;
load('Bout_hom_eta.mat');

V_SD0=grdSp_SD0^3;
V=grdSp^3;

Psrc=75e-3; %W

medNoiPow=0.05;
srcNoiPow=(1-medNoiPow)/2;
detNoiPow=srcNoiPow;
Nnoi=1e6;

%% Calc W
W0=PHI_SD0(:, :, :, 2).*PHI_SD0(:, :, :, 1)*V_SD0;

W1A=PHI(:, :, :, 2).*PHI(:, :, :, 1)*V;
W1B=PHI(:, :, :, 2).*PHI(:, :, :, 4)*V;
W2A=PHI(:, :, :, 3).*PHI(:, :, :, 1)*V;
W2B=PHI(:, :, :, 3).*PHI(:, :, :, 4)*V;

%% Make Noise
C1=1+randn(Nnoi, 1)*noi2bck*sqrt(srcNoiPow);
C2=1+randn(Nnoi, 1)*noi2bck*sqrt(srcNoiPow);
CA=1+randn(Nnoi, 1)*noi2bck*sqrt(detNoiPow);
CB=1+randn(Nnoi, 1)*noi2bck*sqrt(detNoiPow);

C1A=1+randn(Nnoi, 1)*noi2bck*sqrt(medNoiPow);
C1B=1+randn(Nnoi, 1)*noi2bck*sqrt(medNoiPow);
C2A=1+randn(Nnoi, 1)*noi2bck*sqrt(medNoiPow);
C2B=1+randn(Nnoi, 1)*noi2bck*sqrt(medNoiPow);

noi0=1+randn(Nnoi, 1)*noi2bck;
noi1A=C1.*CA.*C1A;
noi1B=C1.*CB.*C1B;
noi2A=C2.*CA.*C2A;
noi2B=C2.*CB.*C2B;

%% R [W/mm^2]
R0bck=Psrc*sum(afDist(:).*W0(:))*etaMuaBck;
R1Abck=Psrc*sum(afDist(:).*W1A(:))*etaMuaBck;
R1Bbck=Psrc*sum(afDist(:).*W1B(:))*etaMuaBck;
R2Abck=Psrc*sum(afDist(:).*W2A(:))*etaMuaBck;
R2Bbck=Psrc*sum(afDist(:).*W2B(:))*etaMuaBck;

R0bck_het=Psrc*sum(afDist_het(:).*W0(:))*etaMuaBck_het;
R1Abck_het=Psrc*sum(afDist_het(:).*W1A(:))*etaMuaBck_het;
R1Bbck_het=Psrc*sum(afDist_het(:).*W1B(:))*etaMuaBck_het;
R2Abck_het=Psrc*sum(afDist_het(:).*W2A(:))*etaMuaBck_het;
R2Bbck_het=Psrc*sum(afDist_het(:).*W2B(:))*etaMuaBck_het;

R0pck=Psrc*W0*etaMuaPck;
R1Apck=Psrc*W1A*etaMuaPck;
R1Bpck=Psrc*W1B*etaMuaPck;
R2Apck=Psrc*W2A*etaMuaPck;
R2Bpck=Psrc*W2B*etaMuaPck;

R0sig=R0pck+R0bck;
R1Asig=R1Apck+R1Abck;
R1Bsig=R1Bpck+R1Bbck;
R2Asig=R2Apck+R2Abck;
R2Bsig=R2Bpck+R2Bbck;

R0sig_het=R0pck+R0bck_het;
R1Asig_het=R1Apck+R1Abck_het;
R1Bsig_het=R1Bpck+R1Bbck_het;
R2Asig_het=R2Apck+R2Abck_het;
R2Bsig_het=R2Bpck+R2Bbck_het;

%% Noise on bck
R0bck_noi=R0bck*noi0;
R1Abck_noi=R1Abck*noi1A;
R1Bbck_noi=R1Bbck*noi1B;
R2Abck_noi=R2Abck*noi2A;
R2Bbck_noi=R2Bbck*noi2B;

R0bck_het_noi=R0bck_het*noi0;
R1Abck_het_noi=R1Abck_het*noi1A;
R1Bbck_het_noi=R1Bbck_het*noi1B;
R2Abck_het_noi=R2Abck_het*noi2A;
R2Bbck_het_noi=R2Bbck_het*noi2B;

%% Sim SD
SD0snr=(R0sig-R0bck)/std(R0bck_noi);
SD1Asnr=(R1Asig-R1Abck)/std(R2Abck_noi);
SD2Asnr=(R2Asig-R2Abck)/std(R2Abck_noi);

SD0snr_het=(R0sig_het-R0bck_het)./std(R0bck_het_noi);
SD1Asnr_het=(R1Asig_het-R1Abck_het)./std(R1Abck_het_noi);
SD2Asnr_het=(R2Asig_het-R2Abck_het)./std(R2Abck_het_noi);

%% Sim DR
DRA12Bsnr=(sqrt((R2Asig.*R1Bsig)./(R1Asig.*R2Bsig))-...
    sqrt((R2Abck.*R1Bbck)./(R1Abck.*R2Bbck)))/...
    std(sqrt((R2Abck_noi.*R1Bbck_noi)./(R1Abck_noi.*R2Bbck_noi)));

DRA12Bsnr_het=(sqrt((R2Asig_het.*R1Bsig_het)./(R1Asig_het.*R2Bsig_het))-...
    sqrt((R2Abck_het.*R1Bbck_het)./(R1Abck_het.*R2Bbck_het)))/...
    std(sqrt((R2Abck_het_noi.*R1Bbck_het_noi)./ ...
    (R1Abck_het_noi.*R2Bbck_het_noi)));

%% Plot Map
cols=turbo(1000);
cols(1, :)=0;
cols(end, :)=1;
cl=[-1, 1]*11;
cTicks=[cl(1):2:-1, 1:2:cl(2)];

xTicks=-5:2.5:5;
yTicks=0:1:4;

yl=[-1, 4.2];

cLin=linspace(cl(1), cl(2), size(cols, 1));
% cols(abs(cLin)<5, :)=0.5;
cols(abs(cLin)<1, :)=0.5;

[~, yInd_SD0]=min(abs(y_SD0-0));
[~, yInd]=min(abs(y-0));

h=figure(10); clf;
h.Name='SNRmaps';
subaxis(2, 2, 1, 'mt', 0.11, 'mb', 0.09, 'ml', 0.05, 'mr', 0.15,...
    'sh', 0.03, 'sv', 0.06); colormap(cols);
tmp=SD0snr_het;
imagesc(x, z, squeeze(tmp(:, yInd_SD0, :)).'); hold on;
caxis(cl);
contour(x, z, squeeze(tmp(:, yInd_SD0, :)).', cTicks,...
    '-k'); 
plot([0, 0, -0.25, 0, 0.25]-0.05, [0, -1, -0.75, -1, -0.75], '-b');
plot([0, 0, -0.25, 0, 0.25]+0.05, [-1, 0, -0.25, 0, -0.25], '-r');
hold off;
shading flat;
set(gca, 'YDir', 'reverse');
axis equal tight;
ylim(yl);
xlim([-5, 5]);
ylabel('$z$ (mm)',...
    'Interpreter', 'latex');
set(gca, 'XTickLabel', {});
set(gca, 'XTick', xTicks);
set(gca, 'YTick', yTicks);
title('\textbf{(a)} Single-Distance (SD), $\rho=0$ mm',...
    'Interpreter', 'latex');

subaxis(2, 2, 2); colormap(cols);
tmp=SD1Asnr_het;
imagesc(x, z, squeeze(tmp(:, yInd, :)).'); hold on;
caxis(cl);
contour(x, z, squeeze(tmp(:, yInd, :)).', cTicks,...
    '-k'); 
plot([0, 0, -0.25, 0, 0.25]-3.5, [0, -1, -0.75, -1, -0.75], '-b');
plot([0, 0, -0.25, 0, 0.25]-0.5, [-1, 0, -0.25, 0, -0.25], '-r');
hold off;
shading flat;
set(gca, 'YDir', 'reverse');
axis equal tight;
ylim(yl);
xlim([-5, 5]);
set(gca, 'YTickLabel', {});
set(gca, 'XTickLabel', {});
set(gca, 'XTick', xTicks);
set(gca, 'YTick', yTicks);
title('\textbf{(b)} Single-Distance (SD), $\rho=3$ mm',...
    'Interpreter', 'latex');

subaxis(2, 2, 3); colormap(cols);
tmp=SD2Asnr_het;
imagesc(x, z, squeeze(tmp(:, yInd, :)).'); hold on;
caxis(cl);
contour(x, z, squeeze(tmp(:, yInd, :)).', cTicks,...
    '-k'); 
plot([0, 0, -0.25, 0, 0.25]-3.5, [0, -1, -0.75, -1, -0.75], '-b');
plot([0, 0, -0.25, 0, 0.25]+0.5, [-1, 0, -0.25, 0, -0.25], '-r');
hold off;
shading flat;
set(gca, 'YDir', 'reverse');
axis equal tight;
ylim(yl);
xlim([-5, 5]);
ylabel('$z$ (mm)');
xlabel('$x$ (mm)');
set(gca, 'XTick', xTicks);
set(gca, 'YTick', yTicks);
title('\textbf{(c)} Single-Distance (SD), $\rho=4$ mm',...
    'Interpreter', 'latex');

subaxis(2, 2, 4); colormap(cols);
tmp=DRA12Bsnr_het;
imagesc(x, z, squeeze(tmp(:, yInd, :)).'); hold on;
pos=get(gca, "Position");
cb=colorbar;
set(gca, "Position", pos);
caxis(cl);
cb.Ticks=cTicks;
contour(x, z, squeeze(tmp(:, yInd, :)).', cTicks,...
    '-k');
plot([0, 0, -0.25, 0, 0.25]+3.5, [0, -1, -0.75, -1, -0.75], '-b');
plot([0, 0, -0.25, 0, 0.25]-3.5, [0, -1, -0.75, -1, -0.75], '-b');
plot([0, 0, -0.25, 0, 0.25]-0.5, [-1, 0, -0.25, 0, -0.25], '-r');
plot([0, 0, -0.25, 0, 0.25]+0.5, [-1, 0, -0.25, 0, -0.25], '-r');
shading flat;
set(gca, 'YDir', 'reverse');
axis equal tight;
ylim(yl);
xlim([-5, 5]);
ylabel(cb, sprintf('Signal-to-Noise Ratio (SNR)\nfrom Fluorescent Probe in Voxel'),...
    'Interpreter', 'latex');
xlabel('$x$ (mm)',...
    'Interpreter', 'latex');
set(gca, 'YTickLabel', {});
set(gca, 'XTick', xTicks);
set(gca, 'YTick', yTicks);
title('\textbf{(d)} Dual-Ratio (DR), $\rho=[3, 4]$ mm',...
    'Interpreter', 'latex');
cb.Units='centimeters';
cb.Position=cb.Position+[0, 0, 0, 3.85];

sgtitle('\textbf{Top-Weighted Autofluorescence}');

%% Plot Hom Map
cols=turbo(1000);
cols(1, :)=0;
cols(end, :)=1;
cl=[-1, 1]*11;
cTicks=[cl(1):2:-1, 1:2:cl(2)];

xTicks=-5:2.5:5;
yTicks=0:1:4;

yl=[-1, 4.2];

cLin=linspace(cl(1), cl(2), size(cols, 1));
% cols(abs(cLin)<5, :)=0.5;
cols(abs(cLin)<1, :)=0.5;

[~, yInd_SD0]=min(abs(y_SD0-0));
[~, yInd]=min(abs(y-0));

h=figure(11); clf;
h.Name='SNRmaps';
subaxis(2, 2, 1, 'mt', 0.11, 'mb', 0.09, 'ml', 0.05, 'mr', 0.15,...
    'sh', 0.03, 'sv', 0.06); colormap(cols);
tmp=SD0snr;
imagesc(x, z, squeeze(tmp(:, yInd_SD0, :)).'); hold on;
caxis(cl);
contour(x, z, squeeze(tmp(:, yInd_SD0, :)).', cTicks,...
    '-k'); 
plot([0, 0, -0.25, 0, 0.25]-0.05, [0, -1, -0.75, -1, -0.75], '-b');
plot([0, 0, -0.25, 0, 0.25]+0.05, [-1, 0, -0.25, 0, -0.25], '-r');
hold off;
shading flat;
set(gca, 'YDir', 'reverse');
axis equal tight;
ylim(yl);
xlim([-5, 5]);
ylabel('$z$ (mm)',...
    'Interpreter', 'latex');
set(gca, 'XTickLabel', {});
set(gca, 'XTick', xTicks);
set(gca, 'YTick', yTicks);
title('\textbf{(e)} Single-Distance (SD), $\rho=0$ mm',...
    'Interpreter', 'latex');

subaxis(2, 2, 2); colormap(cols);
tmp=SD1Asnr;
imagesc(x, z, squeeze(tmp(:, yInd, :)).'); hold on;
caxis(cl);
contour(x, z, squeeze(tmp(:, yInd, :)).', cTicks,...
    '-k'); 
plot([0, 0, -0.25, 0, 0.25]-3.5, [0, -1, -0.75, -1, -0.75], '-b');
plot([0, 0, -0.25, 0, 0.25]-0.5, [-1, 0, -0.25, 0, -0.25], '-r');
hold off;
shading flat;
set(gca, 'YDir', 'reverse');
axis equal tight;
ylim(yl);
xlim([-5, 5]);
set(gca, 'YTickLabel', {});
set(gca, 'XTickLabel', {});
set(gca, 'XTick', xTicks);
set(gca, 'YTick', yTicks);
title('\textbf{(f)} Single-Distance (SD), $\rho=3$ mm',...
    'Interpreter', 'latex');

subaxis(2, 2, 3); colormap(cols);
tmp=SD2Asnr;
imagesc(x, z, squeeze(tmp(:, yInd, :)).'); hold on;
caxis(cl);
contour(x, z, squeeze(tmp(:, yInd, :)).', cTicks,...
    '-k'); 
plot([0, 0, -0.25, 0, 0.25]-3.5, [0, -1, -0.75, -1, -0.75], '-b');
plot([0, 0, -0.25, 0, 0.25]+0.5, [-1, 0, -0.25, 0, -0.25], '-r');
hold off;
shading flat;
set(gca, 'YDir', 'reverse');
axis equal tight;
ylim(yl);
xlim([-5, 5]);
ylabel('$z$ (mm)');
xlabel('$x$ (mm)');
set(gca, 'XTick', xTicks);
set(gca, 'YTick', yTicks);
title('\textbf{(g)} Single-Distance (SD), $\rho=4$ mm',...
    'Interpreter', 'latex');

subaxis(2, 2, 4); colormap(cols);
tmp=DRA12Bsnr;
imagesc(x, z, squeeze(tmp(:, yInd, :)).'); hold on;
pos=get(gca, "Position");
cb=colorbar;
set(gca, "Position", pos);
caxis(cl);
cb.Ticks=cTicks;
contour(x, z, squeeze(tmp(:, yInd, :)).', cTicks,...
    '-k');
plot([0, 0, -0.25, 0, 0.25]+3.5, [0, -1, -0.75, -1, -0.75], '-b');
plot([0, 0, -0.25, 0, 0.25]-3.5, [0, -1, -0.75, -1, -0.75], '-b');
plot([0, 0, -0.25, 0, 0.25]-0.5, [-1, 0, -0.25, 0, -0.25], '-r');
plot([0, 0, -0.25, 0, 0.25]+0.5, [-1, 0, -0.25, 0, -0.25], '-r');
shading flat;
set(gca, 'YDir', 'reverse');
axis equal tight;
ylim(yl);
xlim([-5, 5]);
ylabel(cb, sprintf('Signal-to-Noise Ratio (SNR)\nfrom Fluorescent Probe in Voxel'),...
    'Interpreter', 'latex');
xlabel('$x$ (mm)',...
    'Interpreter', 'latex');
set(gca, 'YTickLabel', {});
set(gca, 'XTick', xTicks);
set(gca, 'YTick', yTicks);
title('\textbf{(h)} Dual-Ratio (DR), $\rho=[3, 4]$ mm',...
    'Interpreter', 'latex');
cb.Units='centimeters';
cb.Position=cb.Position+[0, 0, 0, 3.85];

sgtitle('\textbf{Homogeneous Autofluorescence}');

%% Plot Sig
z_slice=0.5:0.5:4;
[~, zInds_SD0]=min(abs(z_SD0-z_slice));
[~, zInds]=min(abs(z-z_slice));

Nslice=length(z_slice);

cols=turbo(Nslice);

yl=[-10, 15; -1, 4.5];
xl=[-5, 5];

for j=1:2
    h=figure(100+j); clf;
    h.Name='SIGsim';
    subaxis(2, 2, 1,...
        'sv', 0.07, 'sh', 0.03,...
        'mr', 0.03, 'ml', 0.12, 'mt', 0.08);

    area(x_SD0, ones(size(x_SD0)), ...
        'FaceColor', [0.5, 0.5, 0.5], 'FaceAlpha', 0.25, ...
        'EdgeColor', 'none'); hold on;
    area(x_SD0, -ones(size(x_SD0)), ...
        'FaceColor', [0.5, 0.5, 0.5], 'FaceAlpha', 0.25, ...
        'EdgeColor', 'none');
    labs={};
    ax=[];
    for i=((j-1)*Nslice/2+1):...
            (((j-1)*Nslice/2+1)+Nslice/2-1)
        tmp=SD0snr_het(:, yInd_SD0, zInds(i));
        axTmp=plot(x_SD0, tmp, '-',...
            'Color', cols(i, :));
        ax=[ax, axTmp];

        tmp=SD0snr(:, yInd_SD0, zInds(i));
        plot(x_SD0, tmp, '--',...
            'color', cols(i, :));

        labs=[labs, {sprintf('$z=%.1f$ mm', z_slice(i))}];
    end; hold off;
    xlim(xl);
    ylim(yl(j, :));
    if j==2
        set(gca, 'YTick', -1:4);
    end
    if j==1
        legend(ax, labs, 'Location', 'south',...
            'NumColumns', 2);
    else
        legend(ax, labs, 'Location', 'north',...
            'NumColumns', 2);
    end
    set(gca, 'XTickLabel', {});
    ylabel(sprintf(...
        'Signal-to-Noise Ratio (SNR)\nfrom Fluorescent Probe'),...
        'Interpreter', 'latex');
    if j==1
        title('\textbf{(a)} Single-Distance (SD), $\rho=0$ mm',...
            'Interpreter', 'latex');
    else
        title('\textbf{(e)} Single-Distance (SD), $\rho=0$ mm',...
            'Interpreter', 'latex');
    end


    subaxis(2, 2, 2);
    area(x, ones(size(x)), ...
        'FaceColor', [0.5, 0.5, 0.5], 'FaceAlpha', 0.25, ...
        'EdgeColor', 'none'); hold on;
    area(x, -ones(size(x)), ...
        'FaceColor', [0.5, 0.5, 0.5], 'FaceAlpha', 0.25, ...
        'EdgeColor', 'none');
    labs={};
    ax=[];
    for i=((j-1)*Nslice/2+1):...
            (((j-1)*Nslice/2+1)+Nslice/2-1)
        tmp1=SD1Asnr_het(:, yInd, zInds(i));
        tmp2=SD1Asnr(:, yInd, zInds(i));
        plot(x, tmp1, '-',...
            'Color', cols(i, :)); hold on;
        plot(x, tmp2, '--',...
            'color', cols(i, :));
    end
    ax(1)=plot(NaN, NaN, '-k');
    ax(2)=plot(NaN, NaN, '--k');
    hold off;
    if j==1
        legend(ax, ...
            sprintf('Top-Weighted\nAutofluorescence'), ...
            sprintf('Homogeneous\nAutofluorescence'), 'Location', 'south',...
            'NumColumns', 2);
    else
        legend(ax, ...
            sprintf('Top-Weighted\nAutofluorescence'), ...
            sprintf('Homogeneous\nAutofluorescence'), 'Location', 'north',...
            'NumColumns', 2);
    end
    xlim(xl);
    ylim(yl(j, :));
    if j==2
        set(gca, 'YTick', -1:4);
    end
    set(gca, 'XTickLabel', {});
    set(gca, 'YTickLabel', {});
    if j==1
        title('\textbf{(b)} Single-Distance (SD), $\rho=3$ mm',...
            'Interpreter', 'latex');
    else
        title('\textbf{(f)} Single-Distance (SD), $\rho=3$ mm',...
            'Interpreter', 'latex');
    end


    subaxis(2, 2, 3);
    area(x, ones(size(x)), ...
        'FaceColor', [0.5, 0.5, 0.5], 'FaceAlpha', 0.25, ...
        'EdgeColor', 'none'); hold on;
    area(x, -ones(size(x)), ...
        'FaceColor', [0.5, 0.5, 0.5], 'FaceAlpha', 0.25, ...
        'EdgeColor', 'none');
    labs={};
    ax=[];
    for i=((j-1)*Nslice/2+1):...
            (((j-1)*Nslice/2+1)+Nslice/2-1)
        tmp1=SD2Asnr_het(:, yInd, zInds(i));
        tmp2=SD2Asnr(:, yInd, zInds(i));
        plot(x, tmp1, '-',...
            'Color', cols(i, :)); hold on;
        plot(x, tmp2, '--',...
            'color', cols(i, :));
    end; hold off;
    xlim(xl);
    ylim(yl(j, :));
    if j==2
        set(gca, 'YTick', -1:4);
    end
    ylabel(sprintf(...
        'Signal-to-Noise Ratio (SNR)\nfrom Fluorescent Probe'),...
        'Interpreter', 'latex');
    xlabel('$x$ (mm)');
    if j==1
        title('\textbf{(c)} Single-Distance (SD), $\rho=4$ mm',...
            'Interpreter', 'latex');
    else
        title('\textbf{(g)} Single-Distance (SD), $\rho=4$ mm',...
            'Interpreter', 'latex');
    end

    subaxis(2, 2, 4);
    area(x, ones(size(x)), ...
        'FaceColor', [0.5, 0.5, 0.5], 'FaceAlpha', 0.25, ...
        'EdgeColor', 'none'); hold on;
    area(x, -ones(size(x)), ...
        'FaceColor', [0.5, 0.5, 0.5], 'FaceAlpha', 0.25, ...
        'EdgeColor', 'none');
    labs={};
    ax=[];
    for i=((j-1)*Nslice/2+1):...
            (((j-1)*Nslice/2+1)+Nslice/2-1)
        tmp1=DRA12Bsnr_het(:, yInd, zInds(i));
        tmp2=DRA12Bsnr(:, yInd, zInds(i));
        plot(x, tmp1, '-',...
            'Color', cols(i, :)); hold on;
        plot(x, tmp2, '--',...
            'color', cols(i, :));
    end; hold off;
    xlim(xl);
    ylim(yl(j, :));
    if j==2
        set(gca, 'YTick', -1:4);
    end
    set(gca, 'YTickLabel', {});
    xlabel('$x$ (mm)');
    if j==1
        title('\textbf{(d)} Dual-Ratio (DR), $\rho=[3, 4]$ mm',...
            'Interpreter', 'latex');
    else
        title('\textbf{(h)} Dual-Ratio (DR), $\rho=[3, 4]$ mm',...
            'Interpreter', 'latex');
    end
end

%% Plot Sig
cols=lines(4);
ax=[];

yl=[0, 10];
xl=[0, 4.2];

h=figure(200); clf;
h.Name='SNRvZ';

area(z, ones(size(z)), ...
    'FaceColor', [0.5, 0.5, 0.5], 'FaceAlpha', 0.25, ...
    'EdgeColor', 'none'); hold on;

[~, xInd]=min(abs(x_SD0-0));
tmp=squeeze(SD0snr_het(xInd, yInd_SD0, :));
ax(1)=plot(z_SD0, tmp, ...
    'color', cols(1, :)); hold on; 
tmp=squeeze(SD0snr(xInd, yInd_SD0, :));
plot(z_SD0, tmp, '--', ...
    'color', cols(1, :)); 

[~, xInd]=min(abs(x--2));
tmp=squeeze(SD1Asnr_het(xInd, yInd, :));
ax(2)=plot(z, tmp, ...
    'color', cols(2, :));
tmp=squeeze(SD1Asnr(xInd, yInd, :));
plot(z, tmp, '--', ...
    'color', cols(2, :));

[~, xInd]=min(abs(x--1.5));
tmp=squeeze(SD2Asnr_het(xInd, yInd, :));
ax(3)=plot(z, tmp, ...
    'color', cols(3, :));
tmp=squeeze(SD2Asnr(xInd, yInd, :));
plot(z, tmp, '--', ...
    'color', cols(3, :));

[~, xInd]=min(abs(x-0));
tmp=squeeze(DRA12Bsnr_het(xInd, yInd, :));
ax(4)=plot(z, tmp, ...
    'color', cols(4, :));
tmp=squeeze(DRA12Bsnr(xInd, yInd, :));
plot(z, tmp, '--', ...
    'color', cols(4, :));

ax(5)=plot(NaN, NaN, '-k');
ax(6)=plot(NaN, NaN, '--k');
hold off;

xlabel('$z$ (mm)');
ylabel(sprintf(...
    'Signal-to-Noise Ratio (SNR)\nfrom Fluorescent Probe'),...
    'Interpreter', 'latex');

ylim(yl);
xlim(xl);

legend(ax,...
    sprintf('Single-Distance (SD) $0$ mm\n$x=%.1f$ mm', 0), ...
    sprintf('Single-Distance (SD) $3$ mm\n$x=%.1f$ mm', -2), ...
    sprintf('Single-Distance (SD) $4$ mm\n$x=%.1f$ mm', -1.5), ...
    sprintf('Dual-Ratio (DR) $[3,4]$ mm\n$x=%.1f$ mm', 0), ...
    sprintf('Top-Weighted\nAutofluorescence'), ...
    sprintf('Homogeneous\nAutofluorescence'));