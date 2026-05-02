%% Setup
home; clear;

load('MCout_SD0.mat', 'x', 'y', 'z', 'ZZ', 'PHI', 'grdSp');
PHI_SD0=PHI;
x_SD0=x; y_SD0=y; z_SD0=z;
grdSp_SD0=grdSp;
ZZ_SD0=ZZ;
load('MCout.mat', 'x', 'y', 'z', 'ZZ', 'PHI', 'grdSp');

load('Bout_PckNoi.mat');
load('Bout_het_eta.mat');
etaMuaBck_het=etaMuaBck;
afDist_het=afDist;
load('Bout_hom_eta.mat');

V_SD0=grdSp_SD0^3;
V=grdSp^3;

Psrc=75e-3; %W

%% Calc W
W0=PHI_SD0(:, :, :, 2).*PHI_SD0(:, :, :, 1)*V_SD0;

W1A=PHI(:, :, :, 2).*PHI(:, :, :, 1)*V;
W1B=PHI(:, :, :, 2).*PHI(:, :, :, 4)*V;
W2A=PHI(:, :, :, 3).*PHI(:, :, :, 1)*V;
W2B=PHI(:, :, :, 3).*PHI(:, :, :, 4)*V;

%% Loop
medNoiPow_all=linspace(0, 0.2, 100);

SD0_snr1z=NaN(size(medNoiPow_all));
SD1A_snr1z=NaN(size(medNoiPow_all));
SD2A_snr1z=NaN(size(medNoiPow_all));
DRA12A_snr1z=NaN(size(medNoiPow_all));
SD0_het_snr1z=NaN(size(medNoiPow_all));
SD1A_het_snr1z=NaN(size(medNoiPow_all));
SD2A_het_snr1z=NaN(size(medNoiPow_all));
DRA12A_het_snr1z=NaN(size(medNoiPow_all));
for i=1:length(medNoiPow_all)
    medNoiPow=medNoiPow_all(i);
    srcNoiPow=(1-medNoiPow)/2;
    detNoiPow=srcNoiPow;
    Nnoi=1e6;
    rng(100);

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
    
    % Find Deepest
    SD0_het_snr1z(i)=max(ZZ_SD0(SD0snr_het>1));
    SD1A_het_snr1z(i)=max(ZZ(SD1Asnr_het>1));
    SD2A_het_snr1z(i)=max(ZZ(SD2Asnr_het>1));
    DRA12A_het_snr1z(i)=max(ZZ(abs(DRA12Bsnr_het)>1));

    SD0_snr1z(i)=max(ZZ_SD0(SD0snr>1));
    SD1A_snr1z(i)=max(ZZ(SD1Asnr>1));
    SD2A_snr1z(i)=max(ZZ(SD2Asnr>1));
    DRA12A_snr1z(i)=max(ZZ(abs(DRA12Bsnr)>1));
end

%% Plot
cols=lines(4);

ax=[];

h=figure(2000); clf;
h.Name='MaxZ';
ax(1)=plot(medNoiPow_all, SD0_het_snr1z, '-', 'color', cols(1, :));
hold on;
ax(2)=plot(medNoiPow_all, SD1A_het_snr1z, '-', 'color', cols(2, :));
ax(3)=plot(medNoiPow_all, SD2A_het_snr1z, '-', 'color', cols(3, :));
ax(4)=plot(medNoiPow_all, DRA12A_het_snr1z, '-', 'color', cols(4, :));

plot(medNoiPow_all, SD0_snr1z, '--', 'color', cols(1, :));
plot(medNoiPow_all, SD1A_snr1z, '--', 'color', cols(2, :));
plot(medNoiPow_all, SD2A_snr1z, '--', 'color', cols(3, :));
plot(medNoiPow_all, DRA12A_snr1z, '--', 'color', cols(4, :));

ax(5)=plot(NaN, NaN, '-k');
ax(6)=plot(NaN, NaN, '--k');
hold off;

xlabel('Fraction of Non-Cancelable Noise ($p_{NC}$)');
ylabel(sprintf(['Deepest $z$ with Singnal-to-Noise Greater than one:\n' ...
    'max[$z$[SNR$>$1]] (mm)']));

ylim([0, 5]);

leg=legend(ax,...
    sprintf('Single-Distance (SD) $0$ mm'), ...
    sprintf('Single-Distance (SD) $3$ mm'), ...
    sprintf('Single-Distance (SD) $4$ mm'), ...
    sprintf('Dual-Ratio (DR) $[3,4]$ mm'), ...
    sprintf('Top-Weighted\nAutofluorescence'), ...
    sprintf('Homogeneous\nAutofluorescence'), ...
    'NumColumns', 2);
leg.Position=leg.Position+[0.1, 0.1, 0, 0];
