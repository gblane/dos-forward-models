function T = Tslab(rho, s, optProp, m_max)
% T = Tslab(rho, s, optProp, m_max)
%   Continuous-wave transmittance Green's function for a homogeneous slab with
%   extrapolated-boundary conditions, evaluated by the method of images.
%   Equation 4.35 in F. Martelli, S. Del Bianco, A. Ismaelli, and G. Zaccanti,
%   "Light Propagation through Biological Tissue and Other Diffusive Media,"
%   SPIE, 2010, doi:10.1117/3.824746.
%
%   Inputs:
%       rho     - Source-detector distance [mm]
%       s       - Slab thickness [mm]
%       optProp - Struct of optical properties with fields:
%                   .nin  - Index of refraction inside the slab
%                   .nout - Index of refraction outside the slab
%                   .musp - Reduced scattering coefficient [1/mm]
%                   .mua  - Absorption coefficient [1/mm]
%                 (default properties are used, with a warning, if empty)
%       m_max   - Number of image-source pairs each side (default 100)
%   Output:
%       T       - Transmittance Green's function [1/mm^2]
%
%   Requires n2A (this repo, src/analytical) on the path.
    arguments (Input)
        rho (1,1) double; % mm -- Source-detector distance
        s (1,1) double; % mm -- Slab thickness

        optProp struct = [];
        m_max (1,1) double = 100;
    end
    arguments (Output)
        T (1,1) double; % 1/mm^2 -- Transmittance Green's function
    end

    % Set default optical properties if none given
    if isempty(optProp)
        clear optProp;

        optProp.nin = 1.4;
        optProp.nout = 1;
        optProp.musp = 1.1; % 1/mm
        optProp.mua = 0.011; % 1/mm

        warning("Tslab:defaultOptProp", "Default optical properties used");
    end

    mua = optProp.mua; % 1/mm
    D = 1/(3*optProp.musp); % mm

    zs = 1/optProp.musp; % mm
    A = n2A(optProp.nin, optProp.nout);
    ze = 2*A*D; % mm

    m = (-m_max:m_max).';
    z1m = (1-2*m)*s - 4*m*ze - zs; % mm
    z2m = (1-2*m)*s - (4*m-2)*ze + zs; % mm

    T = (1/(4*pi)) * sum( ...
        z1m.*(rho^2+z1m.^2).^(-3/2) .* ...
        (1+(mua*(rho^2+z1m.^2)/D).^(1/2)) .* ...
        exp(-(mua*(rho^2+z1m.^2)/D).^(1/2)) - ...
        z2m.*(rho^2+z2m.^2).^(-3/2) .* ...
        (1+(mua*(rho^2+z2m.^2)/D).^(1/2)) .* ...
        exp(-(mua*(rho^2+z2m.^2)/D).^(1/2)) ...
        ); % 1/mm^2

end
