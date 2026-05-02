function [Y] = twoLayEffHomoOptProp(X, rho, opts)
% [muaHomo, muspHomo] = twoLayEffHomoOptProp(X, rho, opts)
% Giles Blaney Ph.D. Summer 2022
% Inputs:
%   X       - numPropSets X 5 matrix where:
%               X(:, 1): Top layer absorption coefficient (1/mm) 
%               X(:, 2): Bottom layer absorption coefficient (1/mm) 
%               X(:, 3): Top layer reduced scattering coefficient (1/mm) 
%               X(:, 4): Bottom layer reduced scattering coefficient (1/mm)
%               X(:, 5): Top layer thickness (mm)
%   rho     - 1 X numDist vector of source-detector distnaces (mm)
%   opt     - Struct with feilds (Optional [default]):
%               fmod  - Modulation frequency [140.625e6] (Hz)
%               ni    - Index of refraction inside [1.4]
%               no    - Index of refraction outside [1]
%               B     - Radius of 2-layer cylinder [300] (mm)
%               h_end - Number of Bessel function zeros [3000]
% Outputs:
%   Y       - numPropSets X 2 matrix where:
%               Y(:, 1): Effective homogeneous absorption coefficient
%                        (1/mm) 
%               Y(:, 2): Effective homogeneous reduced scattering 
%                        coefficient (1/mm)

    if nargin<=2
        opts.fmod=140.625e6; %Hz
        opts.ni=1.4;
        opts.no=1;
        opts.B=200; %mm
        opts.h_end=3000; 
    end
    
    %% Forward
    R=R_2L_withPreCom(X, rho, opts);
    
    %% Inverse
    recOpts.mueff_tol=1e-4; %1/mm
    recOpts.n_max=10;
    recOpts.omega=2*pi*opts.fmod; %rad/sec
    recOpts.nin=opts.ni;
    recOpts.nout=opts.no;
    recOpts.v=2.99792458e11/recOpts.nin; %mm/sec
    
    Y=NaN(size(R, 1), 2);
    for i=1:size(R, 1)
        [Y(i, 1), Y(i, 2), ~]=...
            DSR2muamuspEB_iterRecov([rho, rho], [R(i, :), R(i, :)],...
            recOpts);
    end

end