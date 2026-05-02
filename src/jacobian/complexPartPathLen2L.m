function [l] = complexPartPathLen2L(rs, r, rd, V, thk, en, optProp, opts)
% Giles Blaney Spring 2020
% [l] = complexPartPathLen2(rs, r, rd, V, thk, en, optProp, opts)
% Inputs:
%   rs      - Source corrdinates. (mm)
%   r       - Center corrdinate of volume. (mm)
%   rd      - Detector corrdinates. (mm)
%   V       - Volume. (mm^3)
%   thk     - Layer thickness. (mm)
%   en      - (OPTIONAL) Bessel function roots.
%   optProp - (OPTIONAL) Struct of optical properties with the following
%             fields:
%                nin  - (default=[1.4, 1.4]) Index of refraction inside. (-)
%                nout - (default=1) Index of refraction outside. (-)
%                musp - (default=[1.20, 0.25] 1/mm) Reduced scattering.
%                       (1/mm)
%                mua  - (default=[0.008, 0.020] 1/mm) Absorption. (1/mm)
%   opts    - (OPTIONAL) Options structure with the following feilds:
%               fmod   - (default=1.40625 Hz) Modulation frequency {Hz}
%               h_end  - (default=2000) Number of Bessel function zeros
%               B      - (default=150 mm) Radius of cylindrical boundary {mm}
% Outputs:
%   l       - Complex partial pathlength. (mm)

    if nargin<=5
        load('zeroOrdBesselRoots.mat');
        
        optProp.nin=[1.4, 1.4];
        optProp.nout=1;
        optProp.musp=[1.20, 0.25]; %1/mm
        optProp.mua=[0.008, 0.020]; %1/mm
        
        opts.fmod=1.40625e8; %Hz
        opts.h_end=2000;
        opts.B=150; %mm
    end
    
    % Place source at origin of cylindrical geometry
    rss=rs-[0, 0, 1/optProp.musp(1)];
    r_cyl=r-rss;
    rd_cyl=rd-rss;
    rs_cyl=rs-rss;
%     rss_cyl=rss-rss;
%     
%     % Shifted pert location for reciprocity
%     rp_cyl=r_cyl-rd_cyl;
    
    PHIrs_r=complexFluence2L(rs_cyl, r_cyl, thk, en, optProp, opts, 'PHIrs_r');
%     PHIrss_rp=complexFluence2L(rss_cyl, rp_cyl, thk, en, optProp, opts,  'PHIrss_rp');
    PHIrd_r=complexFluence2L(rd_cyl, r_cyl, thk, en, optProp, opts, 'PHIrd_r');
    PHIrs_rd=complexFluence2L(rs_cyl, rd_cyl, thk, en, optProp, opts);
    
%     l=(PHIrs_r.*PHIrss_rp.*V)./PHIrs_rd;
    l=(PHIrs_r.*PHIrd_r.*V)./PHIrs_rd;
end