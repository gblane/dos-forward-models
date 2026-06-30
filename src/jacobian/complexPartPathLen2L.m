function [l] = complexPartPathLen2L(rs, r, rd, V, thk, en, optProp, opts)
% complexPartPathLen2L Calculate complex partial pathlength in a two-layer medium.
%
% [l] = complexPartPathLen2L(rs, r, rd, V, thk, en, optProp, opts)
%
% Written by Giles Blaney (Spring 2020; Ph.D. awarded May 2022)
%
% Inputs:
%   rs      - Source coordinates [mm]
%   r       - Center coordinate of volume [mm]
%   rd      - Detector coordinates [mm]
%   V       - Volume [mm^3]
%   thk     - Layer thickness [mm]
%   en      - Bessel function roots [unitless]
%   optProp - Struct of optical properties [struct]
%   opts    - Options structure [struct]
%
% Outputs:
%   l - Complex partial pathlength [mm]

    if nargin<=5
        optProp.nin=[1.4, 1.4];
        optProp.nout=1;
        optProp.musp=[1.20, 0.25]; %1/mm
        optProp.mua=[0.008, 0.020]; %1/mm
        
        opts.fmod=1.40625e8; %Hz
        opts.h_end=2000;
        opts.B=150; %mm
        en=zeroOrdBesselRoots(opts.h_end);
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
