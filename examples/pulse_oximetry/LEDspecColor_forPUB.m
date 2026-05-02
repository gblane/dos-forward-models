%% Giles Blaney, PhD; Spring 2026
clear; home;

%  Nominal RED and IR wavelengths & full width half maxes
lamRED = 660; %nm
lamIR = 880; %nm
LEDfwhm_all = linspace(10,100, 10); %nm

% Melanin volume fractions
Mfrac_all = linspace(0,0.43,10);

% Pulsatile hemoglobin perturbation & arterial oxygen saturation
dTart = 1; %uM
SaO2_all = linspace(0.7,1,100);

% Wavelengths to simulate
lam = (550:1050).'; %nm

% Set if we assume constant epidermis thickness (useMCpathlen=false) or use
% path-lengths from Monte-Carlo (useMCpathlen=true)
useMCpathlen = true;

% Thickness and source detector distance
s = 10; %mm
rho = s; % Detector is directly opposite source

% Number of terms in diffusion theory solution to sum
m_max = 100;

%% ########## Setup #######################################################
% Load results from Monte Carlo based on the model in:
% G. Blaney, J. Frias, F. Tavakoli, A. Sassaroli, and S. Fantini, 
% "Dual-ratio approach to pulse oximetry and the effect of skin tone," 
% JBO, vol. 29, no. S3, p. S33311, Oct. 2024, 
% doi: 10.1117/1.JBO.29.S3.S33311
if useMCpathlen
    MCavg = load('deps/MCout_avg.mat');
end

% Find indexes of the nominal RED and IR wavelengths
[~,iRED] = min(abs(lam-lamRED));
[~,iIR] = min(abs(lam-lamIR));

% Create Laser Diode (LD) emission spectrum
P0red_LD = zeros(size(lam));
P0red_LD(iRED) = 1; %power/nm
P0ir_LD = zeros(size(lam));
P0ir_LD(iIR) = 1; %power/nm

% Generate bulk tissue optical property spectra
[muspTIS, muaTIS, nTIS] = tissueOptProps_func(lam);

% Generate matrix of extinction coefficients for oxy- and deoxy-hemoglobin
E_OD = makeE('OD',lam); %(1/mm)/uM

%% ########## Calculate total path-length in bulk tissue <L> ##############
% Infinitesimal perturbation for numerical derivative 
dmua_numDer = 1e-6; %1/mm

% Initialize optical property structure 
optProp.nin=[];
optProp.nout=1;
optProp.musp=[]; %1/mm
optProp.mua=[]; %1/mm

% Loop through wavelength and calculate <L> for each
L = NaN(size(lam));
for i = 1:length(lam)
    
    % Set refractive index and reduced scattering coefficient
    optProp.nin=nTIS(i);
    optProp.musp=muspTIS(i); %1/mm
    
    % Calculate transmittance Green's function with a symmetric 
    % infinitesimal perturbation
    optProp.mua=muaTIS(i)-dmua_numDer/2; %1/mm
    Tgrn_0 = Tslab(rho, s, optProp, m_max); %1/mm^2
    optProp.mua=muaTIS(i)+dmua_numDer/2; %1/mm
    Tgrn_1 = Tslab(rho, s, optProp, m_max); %1/mm^2
    
    % Find <L> via numerical derivative
    L(i) = -(log(Tgrn_1)-log(Tgrn_0))/dmua_numDer; %mm

end

%% ########## Loop over all arterial oxygen saturations ###################
% Variable initialization (below)

% Intensity spectra emitted from tissue for illumination from the RED or 
% IR LED (wavelength,SaO2,Melanin,FWHM)
IredLED = NaN(length(lam),length(SaO2_all),length(Mfrac_all), ...
    length(LEDfwhm_all));
IirLED = NaN(length(lam),length(SaO2_all),length(Mfrac_all), ...
    length(LEDfwhm_all));

% Ratio-of-ratios for LDs or LEDs (SaO2,Melanin,FWHM)
RoR_LD = NaN(length(SaO2_all),length(Mfrac_all),length(LEDfwhm_all));
RoR_LED = NaN(length(SaO2_all),length(Mfrac_all),length(LEDfwhm_all));

% Recovered arterial oxygen saturations for LDs, LEDs with nominal 
% extinction coefficients, or LEDs with weighted extinction coefficients 
% (SaO2,Melanin,FWHM)
sat_LD_rec = NaN(length(SaO2_all),length(Mfrac_all),length(LEDfwhm_all));
sat_LED_rec = NaN(length(SaO2_all),length(Mfrac_all),length(LEDfwhm_all));
satw_LED_rec = NaN(length(SaO2_all),length(Mfrac_all),length(LEDfwhm_all));

% SaO2 loop
for SaO2ind = 1:length(SaO2_all)
    % Calculate changes in absorption given total-hemoglobin concentration 
    % and arterial oxygen saturation
    dO = dTart*SaO2_all(SaO2ind); %uM
    dD = dTart-dO; %uM
    dmua = E_OD*[dO;dD]; %1/mm
    
    %% ########## Transmittance Green's functions #########################
    % Loop through wavelength and calculate Green's function during 
    % systolic and diastolic phases for each wavelength
    Tgrn_sys = NaN(size(lam));
    Tgrn_dia = NaN(size(lam));
    for i = 1:length(lam)
        
        % Set refractive index and reduced scattering coefficient
        optProp.nin=nTIS(i);
        optProp.musp=muspTIS(i); %1/mm
        
        % Calculate transmittance Green's function with a symmetric 
        % blood volume perturbation
        optProp.mua=muaTIS(i)+dmua(i)/2; %1/mm
        Tgrn_sys(i) = Tslab(rho, s, optProp, m_max); %1/mm^2
        optProp.mua=muaTIS(i)-dmua(i)/2; %1/mm
        Tgrn_dia(i) = Tslab(rho, s, optProp, m_max); %1/mm^2

    end
    
    %% ########## Loop through melanin concentrations #####################
    for melInd = 1:length(Mfrac_all)
        % Absorption of epidermis from melanosomes using Equation 8 in:
        % S. L. Jacques, "Optical properties of biological tissues: a 
        % review," PMB, vol. 58, no. 11, pp. R37–R61, May 2013, 
        % doi: 10.1088/0031-9155/58/11/r37
        muaEPI = Mfrac_all(melInd) * 51.9*(lam/500).^-3.5; %1/mm

        %% Determine length through epidermal melanin filter
        if useMCpathlen
            % Find index of closest melanin concentration in Monte Carlo 
            % data
            [~,iM_MC] = min(abs( MCavg.Mfracs - Mfrac_all(melInd) ));
            
            % Map partial path-lengths from Monte Carlo onto the simulated 
            % wavelengths
            tmpLam = movmean(MCavg.lams, 10); %nm
            tmpL = movmean(MCavg.l(iM_MC,:,1), 10); %mm
            d = interp1(tmpLam, tmpL, lam, ...
                'linear','extrap'); %mm
        else
            % Assume 0.25 mm thick on each side
            d = 0.25*2; %mm
        end
        
        %% ########## Loop through FWHM cases #############################
        for fwhmInd = 1:length(LEDfwhm_all)
            % Generate LED emission spectra 
            P0red_LED = LEDspec_func(lam, lamRED, LEDfwhm_all(fwhmInd));
            P0ir_LED = LEDspec_func(lam, lamIR, LEDfwhm_all(fwhmInd));
        
            %% ########## Model spectra reaching the detector #############
            % Spectra for LD and systolic or diastolic cases (will be the 
            % same for every FWHM loop iteration)
            IredLD_sys = P0red_LD .* exp(-muaEPI.*d) .* Tgrn_sys;
            IirLD_sys = P0ir_LD .* exp(-muaEPI.*d) .* Tgrn_sys;
            IredLD_dia = P0red_LD .* exp(-muaEPI.*d) .* Tgrn_dia;
            IirLD_dia = P0ir_LD .* exp(-muaEPI.*d) .* Tgrn_dia;
            
            % Spectra for LED and systolic or diastolic cases (will depend 
            % on FWHM loop iteration)
            IredLED_sys = P0red_LED .* exp(-muaEPI.*d) .* Tgrn_sys;
            IirLED_sys = P0ir_LED .* exp(-muaEPI.*d) .* Tgrn_sys;
            IredLED_dia = P0red_LED .* exp(-muaEPI.*d) .* Tgrn_dia;
            IirLED_dia = P0ir_LED .* exp(-muaEPI.*d) .* Tgrn_dia;

            % Find average spectra reaching the detectors for the LED case
            IredLED(:,SaO2ind,melInd,fwhmInd) = ...
                (IredLED_sys+IredLED_dia)/2;
            IirLED(:,SaO2ind,melInd,fwhmInd) = ...
                (IirLED_sys+IirLED_dia)/2;
            
            %% ########## Find detected signals ###########################
            % Sum spectra reaching the detectors to find detected signal 
            % for LDs during systolic and diastolic phases  (will be the 
            % same for every FWHM loop iteration)
            IredLD_sys_det = sum(IredLD_sys);
            IirLD_sys_det = sum(IirLD_sys);
            IredLD_dia_det = sum(IredLD_dia);
            IirLD_dia_det = sum(IirLD_dia);
            
            % Sum spectra reaching the detectors to find detected signal 
            % for LEDs during systolic and diastolic phases  (will depend 
            % on FWHM loop iteration)
            IredLED_sys_det = sum(IredLED_sys);
            IirLED_sys_det = sum(IirLED_sys);
            IredLED_dia_det = sum(IredLED_dia);
            IirLED_dia_det = sum(IirLED_dia);
    
            %% ########## Calculate ratio-of-ratios #######################
            RoR_LD(SaO2ind,melInd,fwhmInd) = ...
                ((IredLD_sys_det-IredLD_dia_det)/IredLD_dia_det)/...
                ((IirLD_sys_det-IirLD_dia_det)/IirLD_dia_det);
            RoR_LED(SaO2ind,melInd,fwhmInd) = ...
                ((IredLED_sys_det-IredLED_dia_det)/IredLED_dia_det)/...
                ((IirLED_sys_det-IirLED_dia_det)/IirLED_dia_det);
            
            %% ########## Apply modified Beer Lambert law #################
            % Generate matrix of extinction coefficients for oxy- and 
            % deoxy-hemoglobin only at the nominal wavelengths
            E_mBLL = makeE('OD', [lamRED,lamIR]);
            
            % Calculate recovered absorption changes for LDs
            dmuaRED_LD = -log(IredLD_sys_det/IredLD_dia_det)/L(iRED);
            dmuaIR_LD = -log(IirLD_sys_det/IirLD_dia_det)/L(iIR);
            
            % Convert absorption to concentration and calculate recovered 
            % oxygen saturation for LDs
            dC_LD = E_mBLL\[dmuaRED_LD;dmuaIR_LD];
            dO_LD_rec = dC_LD(1);
            dD_LD_rec = dC_LD(2);
            sat_LD_rec(SaO2ind,melInd,fwhmInd) = ...
                dO_LD_rec / (dO_LD_rec+dD_LD_rec);
            
            % Calculate recovered absorption changes for LEDs
            dmuaRED_LED = -log(IredLED_sys_det/IredLED_dia_det)/L(iRED);
            dmuaIR_LED = -log(IirLED_sys_det/IirLED_dia_det)/L(iIR);
            
            % Convert absorption to concentration and calculate recovered 
            % oxygen saturation for LDs
            dC_LED = E_mBLL\[dmuaRED_LED;dmuaIR_LED];
            dO_LED_rec = dC_LED(1);
            dD_LED_rec = dC_LED(2);
            sat_LED_rec(SaO2ind,melInd,fwhmInd) = ...
                dO_LED_rec / (dO_LED_rec+dD_LED_rec);
            
            %% ########## Implement weighted extinction coefficients ######
            % Generate extinction spectra for oxy- and deoxy-hemoglobin 
            extO = makeE('O',lam);
            extD = makeE('D',lam);

            % Weight extinction coefficients by the spectra reaching the 
            % detectors
            Wred = IredLED(:,SaO2ind,melInd,fwhmInd);
            Wir = IirLED(:,SaO2ind,melInd,fwhmInd);
            extOred = sum(Wred.*extO)/sum(Wred);
            extDred = sum(Wred.*extD)/sum(Wred);
            extOir = sum(Wir.*extO)/sum(Wir);
            extDir = sum(Wir.*extD)/sum(Wir);
            
            % Convert absorption to concentration and calculate recovered 
            % oxygen saturation considering weighted coefficients for LEDs
            Ew_mBLL = [
                extOred, extDred;
                extOir, extDir
                ];
            dCw_LED = Ew_mBLL\[dmuaRED_LED;dmuaIR_LED];
            dOw_LED_rec = dCw_LED(1);
            dDw_LED_rec = dCw_LED(2);
            satw_LED_rec(SaO2ind,melInd,fwhmInd) = ...
                dOw_LED_rec / (dOw_LED_rec+dDw_LED_rec);
        end
    end
end

%% ########## Plot -- Modeled Spectra #####################################
cols = lines(7);
Mind = length(Mfrac_all);

[~,satInd] = min(abs(SaO2_all-0.95));
[~,fwhmInd] = min(abs(LEDfwhm_all-50));



figure(101); clf;

% -------------------------------------------------------------------------
subplot(2,1,1); hold on;
plot(lam, LEDspec_func(lam, lamRED, LEDfwhm_all(fwhmInd)), '-',...
    'Color',cols(7,:));
xlim([-100,100]+lamRED);
ylim([0,1]);
ylabel('$P/P_{max}$');
xlabel('Wavelength (nm)');
title('Red LED');

pos = get(gca,"Position");
legend( ...
    'Source Emission', ...
    'Orientation','horizontal', ...
    'Location','northoutside');
set(gca, "Position",pos);

% -------------------------------------------------------------------------
subplot(2,1,2); hold on;
plot(lam, LEDspec_func(lam, lamRED, LEDfwhm_all(fwhmInd)), '--',...
    'Color',cols(7,:));
tmp = IredLED(:,satInd,Mind,fwhmInd);
plot(lam, tmp/max(tmp), '-',...
    'Color',cols(7,:));
xlim([-100,100]+lamRED);
ylim([0,1]);
ylabel('$I/I_{max}$');
xlabel('Wavelength (nm)');
title('Red LED');

pos = get(gca,"Position");
legend( ...
    'Source Emission', ...
    'Detected Spectrum', ...
    'Orientation','vertical', ...
    'Location','northoutside');
set(gca, "Position",pos);



% -------------------------------------------------------------------------
% Absorption of epidermis from melanosomes using Equation 8 in:
% S. L. Jacques, "Optical properties of biological tissues: a 
% review," PMB, vol. 58, no. 11, pp. R37–R61, May 2013, 
% doi: 10.1088/0031-9155/58/11/r37
muaEPI = Mfrac_all(Mind) * 51.9*(lam/500).^-3.5; %1/mm

[~,iM_MC] = min(abs( MCavg.Mfracs - Mfrac_all(Mind) ));
tmpLam = movmean(MCavg.lams, 10);
tmpL = movmean(MCavg.l(iM_MC,:,1), 10);
d = interp1(tmpLam, tmpL, lam, ...
    'linear','extrap');

figure(102); clf;
yyaxis left;
plot(lam, muaEPI);
ylabel('$\mu_a$ (mm$^{-1}$)');

yyaxis right;
plot(lam, d);
ylabel('$\langle \ell \rangle$ (mm)', ...
    'Interpreter','latex');

xlim(lam([1,end]));
xlabel('Wavelength (nm)');

title('Epidermis Properties');



% -------------------------------------------------------------------------
figure(103); clf;

yyaxis left;
plot(lam, muaTIS);
ylabel('$\mu_a$ (mm$^{-1}$)');

yyaxis right;
plot(lam, L);
ylabel('$\langle L \rangle$ (mm)');

xlim(lam([1,end]));
xlabel('Wavelength (nm)');

title('Bulk Tissue Properties');

%% ########## Plot -- SaO2 vs. RoR ########################################

colsM=flipud(copper(length(Mfrac_all)));

[~,fwhmInd] = min(abs(LEDfwhm_all-50));



figure(200); clf;
% -------------------------------------------------------------------------
ax1 = subplot(2,2,1);
slp_LD = NaN(size(Mfrac_all));
for i = 1:length(Mfrac_all)
    plot(RoR_LD(:,i,fwhmInd), SaO2_all,'-',...
        'Color',colsM(i,:));
    hold on;
    
    p = polyfit(RoR_LD(:,i,fwhmInd),SaO2_all,1);
    slp_LD(i) = p(1);
end
xlabel('Recovered Ratio-of-Ratios');
ylabel('Modeled SaO$_2$');
title('(a) LDs');

% -------------------------------------------------------------------------
ax2 = subplot(2,2,2);
slp_LED = NaN(size(Mfrac_all));
for i = 1:length(Mfrac_all)
    plot(RoR_LED(:,i,fwhmInd), SaO2_all,'-',...
        'Color',colsM(i,:));
    hold on;

    p = polyfit(RoR_LED(:,i,fwhmInd),SaO2_all,1);
    slp_LED(i) = p(1);
end
xlabel('Recovered Ratio-of-Ratios');
ylabel('Modeled SaO$_2$');
title('(b) LEDs');

linkaxes([ax1,ax2], 'xy');
xlim([0.3,1.3]);

% -------------------------------------------------------------------------
ax3 = subplot(2,2,3);
plot(Mfrac_all,slp_LD, '-k'); hold on;
for i = 1:length(Mfrac_all)
    plot(Mfrac_all(i),slp_LD(i), 'o', 'Color',colsM(i,:));
end
xlabel('Melanin Volume Fraction');
ylabel('Slope of (a)');
title('(c) LDs');

% -------------------------------------------------------------------------
ax4 = subplot(2,2,4);
plot(Mfrac_all,slp_LED, '-k'); hold on;
for i = 1:length(Mfrac_all)
    plot(Mfrac_all(i),slp_LED(i), 'o', 'Color',colsM(i,:));
end
xlabel('Melanin Volume Fraction');
ylabel('Slope of (b)');
title('(d) LEDs');

linkaxes([ax3,ax4], 'xy');
ylim([-0.38,-0.3]);

%% ########## Plot -- Application of weighted extinction coefficients #####

colsM = flipud(copper(length(Mfrac_all)));

[~,fwhmInd] = min(abs(LEDfwhm_all-50));



figure(300); clf;
% -------------------------------------------------------------------------
ax1 = subplot(2,2,1);
for i = 1:length(Mfrac_all)
    plot(SaO2_all, sat_LD_rec(:,i,fwhmInd)-SaO2_all.', '-', ...
        'Color',colsM(i,:)); hold on;
end
xlabel('Modeled SaO$_2$');
ylabel('SaO$_{2,recovered}$ $-$ SaO$_{2,modeled}$');
title('(a) LDs with Assumed $\epsilon$s');

% -------------------------------------------------------------------------
ax2 = subplot(2,2,2);
for i = 1:length(Mfrac_all)
    plot(SaO2_all, sat_LED_rec(:,i,fwhmInd)-SaO2_all.', '-', ...
        'Color',colsM(i,:)); hold on;
end
xlabel('Modeled SaO$_2$');
ylabel('SaO$_{2,recovered}$ $-$ SaO$_{2,modeled}$');
title('(b) LEDs with Assumed $\epsilon$s');

% -------------------------------------------------------------------------
subplot(2,2,3);
cb = colorbar;
cb.Limits = Mfrac_all([1,end]);
cb.Ticks = round(Mfrac_all,2);
cb.Label.String = 'Melanin Volume Fraction';
cb.Position = cb.Position + [-0.1, 0, 0.01, 0];
clim(Mfrac_all([1,end]));
colormap(colsM);
axis off;

% -------------------------------------------------------------------------
ax4 = subplot(2,2,4);
for i = 1:length(Mfrac_all)
    plot(SaO2_all, satw_LED_rec(:,i,fwhmInd)-SaO2_all.', '-', ...
        'Color',colsM(i,:)); hold on;
end
xlabel('Modeled SaO$_2$');
ylabel('SaO$_{2,recovered}$ $-$ SaO$_{2,modeled}$');
title('(c) LEDs with Weighted $\epsilon$s');

linkaxes([ax1,ax2,ax4], 'xy');

%% ########## Plot -- Dependance on LED FWHM ##############################

colsM = flipud(copper(length(Mfrac_all)));



figure(400); clf;
% -------------------------------------------------------------------------
ax1 = subplot(1,2,1);
for i = 1:length(Mfrac_all)
    plot([0,LEDfwhm_all], ...
        [squeeze(rms(sat_LD_rec(:,i,1) - SaO2_all.'));...
        squeeze(rms(sat_LED_rec(:,i,:) - SaO2_all.'))], '-', ...
        'Color',colsM(i,:)); hold on;
end
xlabel('Source FWHM (nm)');
ylabel('SaO$_2$ RMS Error');

title('(a) Assumed $\epsilon$s');

% -------------------------------------------------------------------------
ax2 = subplot(1,2,2);
for i = 1:length(Mfrac_all)
    plot([0,LEDfwhm_all], ...
        [squeeze(rms(sat_LD_rec(:,i,1) - SaO2_all.'));...
        squeeze(rms(satw_LED_rec(:,i,:) - SaO2_all.'))], '-', ...
        'Color',colsM(i,:)); hold on;
end
xlabel('Source FWHM (nm)');
ylabel('SaO$_2$ RMS Error');
title('(b) Weighted $\epsilon$s');

pos = get(gca,'Position');
cb = colorbar;
cb.Limits = Mfrac_all([1,end]);
cb.Ticks = round(Mfrac_all,2);
cb.Label.String = 'Melanin Volume Fraction';
clim(Mfrac_all([1,end]));
colormap(colsM);
set(gca,'Position',pos);

linkaxes([ax1,ax2],'xy');



%% ########################################################################
%% ########## Functions ###################################################
%% ########################################################################

% =============== LED Spectra =============================================
function [P0_LED] = LEDspec_func(lam, lamPK, lamFWHM)
    
    % Generate Gaussian and normalize
    P0_LED = normpdf(lam, ...
        lamPK, lamFWHM/(2*sqrt(2*log(2))));
    P0_LED = P0_LED/max(P0_LED);

end

% =============== Tissue Optical Properties ===============================
function [muspTIS, muaTIS, nTIS] = tissueOptProps_func(lam, NVA)
    % Tissue optical properties based on:
    % S. L. Jacques, "Optical properties of biological tissues: a 
    % review," PMB, vol. 58, no. 11, pp. R37–R61, May 2013, 
    % doi: 10.1088/0031-9155/58/11/r37

    % By default use relative tissue volumes from model in:
    % G. Blaney, J. Frias, F. Tavakoli, A. Sassaroli, and S. Fantini, 
    % "Dual-ratio approach to pulse oximetry and the effect of skin tone," 
    % JBO, vol. 29, no. S3, p. S33311, Oct. 2024, 
    % doi: 10.1117/1.JBO.29.S3.S33311
    % Thus V_tisTyp/V_tot,noEpi are:
    % Bone = 0.119
    % Muscle = 0.457
    % Fat = 0.228
    % Dermis = 0.196
    arguments
        lam (:,1) double; %nm
        
        NVA.T (1,1) double = ...
            69.8*0.119 + 117*0.457 + 12.5*0.228 + 4.70*0.196; %uM
        NVA.S (1,1) double = ...
            0.875*0.119 + 0.641*0.457 + 0.760*0.228 + 0.390*0.196;
        NVA.W (1,1) double = ...
            0.318*0.119 + 0.795*0.457 + 0.110*0.228 + 0.650*0.196;
        NVA.L (1,1) double = ...
            0.00*0.119 + 0.00*0.457 + 0.69*0.228 + 0.00*0.196;
    
        NVA.ap (1,1) double = ...
            15.3*0.119 + 13.0*0.457 + 34.2*0.228 + 43.6*0.196; %1/cm
        NVA.fray (1,1) double = ...
            0.022*0.119 + 0.000*0.457 + 0.260*0.228 + 0.410*0.196;
        NVA.bMie (1,1) double = ...
            0.326*0.119 + 0.926*0.457 + 0.567*0.228 + 0.562*0.196;
    end
    
    %% musp
    ap = NVA.ap; % 1/cm
    fRay = NVA.fray;
    bMie = NVA.bMie;
    % Combination of Rayleigh and Mie scattering
    muspTIS = ap*0.1*(fRay*(lam/500).^-4 ...
        + (1-fRay)*(lam/500).^-bMie); %1/mm
    
    %% mua
    T = NVA.T; % uM
    S = NVA.S;
    W = NVA.W;
    L = NVA.L;
    
    % Beer's law
    E_ODWL = makeE('ODWL', lam);
    muaTIS = E_ODWL * [T*S; T*(1-S); W; L];
    
    %% n
    % https://www.engineersedge.com/physics/refraction_for_water__15690.htm
    % at 40 deg C
    n_waterRef=[...
         226.50,  361.05,  404.41,  589.00,  632.80, 1013.98;
        1.39046, 1.34540, 1.34065, 1.33095, 1.32972, 1.32296].';
    n_water=interp1(n_waterRef(:, 1), n_waterRef(:, 2), lam, ...
        'linear', 'extrap');
    
    % Equation 3 in:
    % S. L. Jacques, "Optical properties of biological tissues: a 
    % review," PMB, vol. 58, no. 11, pp. R37–R61, May 2013, 
    % doi: 10.1088/0031-9155/58/11/r37
    n_dry = 1.514;
    nTIS = n_dry - (n_dry - n_water) * W;

end

% =============== Transmittance Green's Function ==========================
function T = Tslab(rho, s, optProp, m_max)
    % Transmitance for a slab with extrapolated boundary
    % Equation 4.35 in:
    % F. Martelli, S. Del Bianco, A. Ismaelli, and G. Zaccanti, Light 
    % Propagation through Biological Tissue and Other Diffusive Media. 
    % 1000 20th Street, Bellingham, WA 98227-0010 USA: SPIE, 2010. 
    % doi: 10.1117/3.824746
    arguments (Input)
        rho (1,1) double; %mm -- Source detector distance 
        s (1,1) double; %mm -- Slab thickness
        
        optProp struct = [];
        m_max (1,1) double = 100;
    end
    arguments (Output)
        T (1,1) double; %1/mm^2 -- Transmittance Green's function
    end
    
    % Set default optical properties if input not given
    if isempty(optProp)
        clear optProp;
        
        optProp.nin=1.4;
        optProp.nout=1;
        optProp.musp=1.1; %1/mm
        optProp.mua=0.011; %1/mm
        
        warning('Default optical properties used');
    end
    
    mua = optProp.mua; %1/mm
    D = 1/(3*optProp.musp); %mm
    
    zs = 1/optProp.musp; %mm
    A = n2A(optProp.nin,optProp.nout);
    ze = 2*A*D; %mm
    
    m = (-m_max:m_max).';
    z1m = (1-2*m)*s - 4*m*ze - zs; %mm
    z2m = (1-2*m)*s - (4*m-2)*ze + zs; %mm
    
    T = (1/(4*pi)) * sum( ...
        z1m.*(rho^2+z1m.^2).^(-3/2) .* ...
        (1+(mua*(rho^2+z1m.^2)/D).^(1/2)) .* ...
        exp(-(mua*(rho^2+z1m.^2)/D).^(1/2)) - ...
        z2m.*(rho^2+z2m.^2).^(-3/2) .* ...
        (1+(mua*(rho^2+z2m.^2)/D).^(1/2)) .* ...
        exp(-(mua*(rho^2+z2m.^2)/D).^(1/2)) ...
        ); %1/mm^2

end

% =============== Reflection Parameter ====================================
function A = n2A(nin, nout)
% Based on:
% R. Aronson, "Boundary conditions for diffusion of light," JOSA A, 
% vol. 12, no. 11, pp. 2532–2539, Nov. 1995, doi: 10/dff63n.

% A = n2A(nin,nout)
%   Inputs:
%       nin  - Index of refraction inside medium
%       nout - Index of refraction outside medium
%   Output:
%       A    - Index of refraction mismatch parameter

    if nargin==0
        dan12=1.4;
    else
        dan12=nin/nout;
    end

    if dan12>1
        A=504.332889-2641.00214*dan12+...
            5923.699064*dan12.^2-7376.355814*dan12^3+...
            5507.53041*dan12^4-2463.357945*dan12^5+...
            610.956547*dan12^6-64.8047*dan12^7;
    elseif dan12<1
        A=3.084635-6.531194*dan12+...
            8.357854*dan12^2-5.082751*dan12^3+1.171382*dan12^4;
    else
        A=1;
    end

end

% =============== Extinction Coefficients ===============
function E = makeE(chroms, lambda)
% E = makeE(chroms, lambda)
% Giles Blaney Spring 2021
% 
% Inputs:   - chroms: String of chromophores to include in E. 
%                     Available chromophores:
%                     - O: Oxyhemoglobin
%                     - D: Deoxyhemoglobin
%                     - W: Water
%                     - L: Lipid
%                     (Default: 'OD')
%                     Spaces are ignored.
%           - lambda: Vectors of wavelengths (nm).
%                     (Default: [830, 690])
% 
% Output:   - E: Extinction coefficient matrix.
%                Units: 1/(mm uM) for O, D, CCOo, and CCOr
%                       1/mm for W, L, and C
%                size(E)=[length(lambda), length(chroms)];
%                Defined as mua=E*C
%                Order of C is defined by order in chroms input
    
% Relevant references in .mat files for each chromophore

    if nargin<=0
        chroms='OD';
    end
    if nargin<=1
        lambda=[830, 690]; %nm
    end
    if size(lambda, 1)==1
        lambda=lambda';
    end
    
    chroms=chroms(~isspace(chroms));
    
    E=[];
    while ~isempty(chroms)
        switch chroms(1)
            case 'O' %Oxy
                blood=load('Bext.mat');
                Oext=interp1(blood.Blambda, blood.Oext, lambda,...
                    'linear', 'extrap');
                Oext(Oext<0)=0;
                E=[E, Oext];
            case 'D' %Deoxy
                blood=load('Bext.mat');
                Dext=interp1(blood.Blambda, blood.Dext, lambda,...
                    'linear', 'extrap');
                Dext(Dext<0)=0;
                E=[E, Dext];
            case 'W' %Water
                water=load('Wext.mat');
                Wext=interp1(water.Wlambda, water.Wmua, lambda,...
                    'linear', 'extrap');
                Wext(Wext<0)=0;
                E=[E, Wext];
            case 'L' %Lipid
                lipid=load('Lext.mat');
                Lext=interp1(lipid.Llambda, lipid.Lmua, lambda,...
                    'linear', 'extrap');
                Lext(Lext<0)=0;
                E=[E, Lext];
            otherwise
                warning(['Unknown chromophore' chroms(1) ', ignored']);
        end
        chroms(1)=[];
    end
end