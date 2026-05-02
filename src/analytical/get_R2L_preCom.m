function [preCom] = get_R2L_preCom(rho, opt)
% [preCom] = get_R2L_preCom(rho, opt)
% Supporting function for R_2L_withPreCom()
% Giles Blaney Ph.D. Summer 2022
% Inputs:
%   rho     - 1 X numDist X 1 vector of source-detector distnaces (mm)
%   opt     - Struct with feilds (Optional [default]):
%               no    - Index of refraction outside [1]
%               ni    - Index of refraction inside [1.4]
%               B     - Radius of 2-layer cylinder [300] (mm)
%               h_end - Number of Bessel function zeros [3000]
% Outputs:
%   preCom  - Struct for use with R_2L_withPreCom() function with feilds:
%               A         - Index of refraction mismatch parameter
%               en_pillar - 1 X 1 X h_end pillar of Bessel function zeros
%               Q         - 1 X numDist X h_end array of Bessel function of
%                           1st kind values
    
    if nargin<=1
        no=1;
        ni=1.4;
        B=300; %mm
        h_end=3000;
    else
        no=opt.no;
        ni=opt.ni;
        B=opt.B; %mm
        h_end=opt.h_end;
    end
    
    load('zeroOrdBesselRoots.mat', 'en');
    en_pillar=reshape(en, 1, 1, length(en));
    
    Q=besselj(...
        0, rho.*en_pillar(1:h_end)/B)./(besselj(1, en_pillar(1:h_end))).^2;
    
    preCom.A=n2A(ni, no);
    preCom.en_pillar=en_pillar(1:h_end);
    preCom.Q=Q;
end