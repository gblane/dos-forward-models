function [l] = temporalVarPartPathLen(rs, r, rd, Vol, optProp, NVA)
% Giles Blaney Ph.D. Spring 2023
% [l] = temporalVarPartPathLen(rs, r, rd, Vol, optProp, NVA)
% Inputs:
%   rs      - Source corrdinates. (mm)
%   r       - Center corrdinate of volume. (mm)
%   rd      - Detector corrdinates. (mm)
%   Vol     - Volume. (mm^3)
%   optProp - (OPTIONAL) Struct of optical properties with the following
%             fields:
%                nin  - (default=1.4) Index of refraction inside. (-)
%                nout - (default=1) Index of refraction outside. (-)
%                musp - (default=1.2 1/mm) Reduced scattering. (1/mm)
%                mua  - (default=0.01 1/mm) Absorption. (1/mm)
%   Name Value Arguments:
%           - 'conv_t' (default=10e3 ps): Time window for convolution. (ps)
%           - 'conv_dt' (default=1 ps): Time step for convolution. (ps)
%           - 'usePar' (default=true): Use parfoor loops.
%           - 'FFTconv' (default=true): Use ifft(fft*fft) as conv.
% Outputs:
%   l       - Partial pathlength of variance. (mm)
    
    arguments
        rs (1,3) double; %mm
        r (:,3) double; %mm
        rd (1,3) double; %mm
        Vol (1,1) double; %mm^3

        optProp struct = [];

        NVA.conv_t (1,1) double = 10e3; %ps
        NVA.conv_dt (1,1) = 1; %ps
        NVA.usePar (1,1) logical = true;
        NVA.FFTconv (1,1) logical = true;
    end
        
    if isempty(optProp)
        clear optProp;
        
        optProp.nin=1.4;
        optProp.nout=1;
        optProp.musp=1.2; %1/mm
        optProp.mua=0.01; %1/mm
        
        warning(['Default optical properties used, this may be inconsistent'...
            ' with the musp used for source depth']);
    end
    
    if size(rs, 1)>1 && size(rd, 1)>1
        error('Can not use multiple sources and multiple detectors');
    end

    NVAstruct=struct2pairs(NVA);
    
    t1Mom=temporalKthMoment(rs, rd, 1, optProp);
    t2Mom=temporalKthMoment(rs, rd, 2, optProp);
    V=t2Mom-t1Mom.^2;

    lt1=temporalKthMomPartPathLen(rs, r, rd, Vol, 1, optProp, NVAstruct{:});
    lt2=temporalKthMomPartPathLen(rs, r, rd, Vol, 2, optProp, NVAstruct{:});
    
    l=(t2Mom.*lt2-2*t1Mom.^2.*lt1)./V;

    l(isnan(l))=0;
end