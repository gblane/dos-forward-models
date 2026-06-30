function [muspTIS, muaTIS, nTIS] = tissueOptProps_func(lam, NVA)
% [muspTIS, muaTIS, nTIS] = tissueOptProps_func(lam, NVA)
%   Bulk soft-tissue optical properties versus wavelength: reduced scattering
%   (Rayleigh + Mie), absorption (oxy/deoxy hemoglobin + water + lipid, via
%   makeE), and refractive index. By default the chromophore amounts and
%   scattering parameters are the relative-tissue-volume model of Blaney et al.
%   (JBO 2024); override any of them through the name-value arguments to model a
%   specific tissue (e.g. a calibrated dermis).
%
%   Inputs:
%       lam      - Wavelength grid [nm], column vector
%     Name-value arguments (all optional; defaults from Blaney et al. 2024):
%       T        - Total hemoglobin concentration [uM]
%       S        - Hemoglobin oxygen saturation [-]
%       W        - Water volume fraction [-]
%       L        - Lipid volume fraction [-]
%       ap       - Reduced-scattering amplitude [1/cm]
%       fray     - Rayleigh fraction [-]
%       bMie     - Mie scattering power [-]
%   Outputs:
%       muspTIS  - Reduced scattering coefficient [1/mm]
%       muaTIS   - Absorption coefficient [1/mm]
%       nTIS     - Refractive index [-]
%
%   Requires makeE (dos-inverse-models) on the path for the chromophore
%   extinction spectra.
%
%   References:
%     S. L. Jacques, "Optical properties of biological tissues: a review," PMB,
%       vol. 58, no. 11, pp. R37-R61, 2013, doi:10.1088/0031-9155/58/11/r37
%     G. Blaney, J. Frias, F. Tavakoli, A. Sassaroli, and S. Fantini, "Dual-ratio
%       approach to pulse oximetry and the effect of skin tone," JBO, vol. 29,
%       no. S3, p. S33311, 2024, doi:10.1117/1.JBO.29.S3.S33311
%   The default chromophore set uses relative tissue volumes V_tisTyp/V_tot,noEpi
%   (Bone 0.119, Muscle 0.457, Fat 0.228, Dermis 0.196).
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

    %% Reduced scattering: combination of Rayleigh and Mie (Jacques 2013)
    ap = NVA.ap; % 1/cm
    fRay = NVA.fray;
    bMie = NVA.bMie;
    muspTIS = ap*0.1*(fRay*(lam/500).^-4 ...
        + (1-fRay)*(lam/500).^-bMie); %1/mm

    %% Absorption: Beer's law over [HbO2; Hb; water; lipid]
    T = NVA.T; % uM
    S = NVA.S;
    W = NVA.W;
    L = NVA.L;
    E_ODWL = makeE('ODWL', lam);
    muaTIS = E_ODWL * [T*S; T*(1-S); W; L]; %1/mm

    %% Refractive index: dry-tissue / water mixture (water n at 40 deg C)
    n_waterRef = [...
         226.50,  361.05,  404.41,  589.00,  632.80, 1013.98;
        1.39046, 1.34540, 1.34065, 1.33095, 1.32972, 1.32296].';
    n_water = interp1(n_waterRef(:, 1), n_waterRef(:, 2), lam, ...
        'linear', 'extrap');
    n_dry = 1.514;                      % Jacques (2013), Eq. 3
    nTIS = n_dry - (n_dry - n_water) * W;

end
