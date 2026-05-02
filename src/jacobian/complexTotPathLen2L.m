function [L, Ll, PHI] = complexTotPathLen2L(rs, rd, thk, en, optProp, opts)
% complexTotPathLen2L Calculate complex total pathlength in a two-layer medium.
%
% [L, Ll, PHI] = complexTotPathLen2L(rs, rd, thk, en, optProp, opts)
%
% Written by Giles Blaney, Ph.D. Spring 2020
%
% Inputs:
%   rs      - Source coordinates [mm]
%   rd      - Detector coordinates [mm]
%   thk     - Layer thickness [mm]
%   en      - Bessel function roots [unitless]
%   optProp - Struct of optical properties [struct]
%   opts    - Options structure [struct]
%
% Outputs:
%   L   - Complex total pathlength [mm]
%   Ll  - Complex pathlength for each layer [mm]
%   PHI - Complex fluence [1/mm^2]

    if nargin<=3
        load('zeroOrdBesselRoots.mat');
        
        optProp.nin=[1.4, 1.4];
        optProp.nout=1;
        optProp.musp=[1.20, 0.25]; %1/mm
        optProp.mua=[0.008, 0.020]; %1/mm
        
        opts.fmod=1.40625e8; %Hz
        opts.h_end=2000;
        opts.B=150; %mm
        opts.Method = 'CompFluence';
    end

    if ~isfield(opts, 'Method')
        opts.Method = 'CompFluence';
    end

    if size(rd, 2)~=3
        rd=rd';
    end
    
    dmua=1e-7; %1/mm
    optProp1=optProp;
    optProp1.mua(1)=optProp.mua(1)+dmua;
    optProp2=optProp;
    optProp2.mua(2)=optProp.mua(2)+dmua;
    
    % Place source at origin of cylindrical geometry
    rss=rs-[0, 0, 1/optProp.musp(1)];
    rd_cyl=rd-rss;
    rs_cyl=rs-rss;
    
    if strcmp('CompFluence',opts.Method)
    PHI=complexFluence2L(rs_cyl, rd_cyl, thk, en, optProp, opts);
    PHI1=complexFluence2L(rs_cyl, rd_cyl, thk, en, optProp1, opts);
    PHI2=complexFluence2L(rs_cyl, rd_cyl, thk, en, optProp2, opts);
    
    elseif strcmp('CompReflectance', opts.Method)
        rho = rs(1) - rd(1); % find a better way to do this 
        X(:,1) = optProp.mua(1);
        X(:,2) = optProp.mua(2);
        X(:,3) = optProp.musp(1);
        X(:,4) = optProp.musp(2); 
        X(:,5) = thk; 
        PHI = R_2L_withPreCom(X, rho, opt); 
        clear X
        X(:,1) = optProp1.mua(1);
        X(:,2) = optProp1.mua(2);
        X(:,3) = optProp1.musp(1);
        X(:,4) = optProp1.musp(2); 
        X(:,5) = thk; 
        PHI1 = R_2L_withPreCom(X, rho, opt); 
        clear X
        X(:,1) = optProp2.mua(1);
        X(:,2) = optProp2.mua(2);
        X(:,3) = optProp2.musp(1);
        X(:,4) = optProp2.musp(2); 
        X(:,5) = thk; 
        PHI2 = R_2L_withPreCom(X, rho, opt); 
    end 
    Ll(:, 1)=-(PHI1-PHI)./(PHI*dmua);
    Ll(:, 2)=-(PHI2-PHI)./(PHI*dmua);
    L=sum(Ll, 2);

end