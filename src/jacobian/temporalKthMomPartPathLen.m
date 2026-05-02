function [l] = temporalKthMomPartPathLen(rs, r, rd, Vol, k, optProp, NVA)
% Giles Blaney Ph.D. Spring 2023
% [l] = temporalKthMomPartPathLen(rs, rd, tns, optProp)
% Inputs:
%   rs      - Source corrdinates. (mm)
%   r       - Center corrdinate of volume. (mm)
%   rd      - Detector corrdinates. (mm)
%   Vol     - Volume. (mm^3)
%   k       - (OPTIONAL; default=1) Moment order.
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
%   l       - Partial pathlength of kth momment of t. (mm)
    
    arguments
        rs (1,3) double; %mm
        r (:,3) double; %mm
        rd (1,3) double; %mm
        Vol (1,1) double; %mm^3

        k (1,1) double = 1; 

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
    
    if NVA.usePar
        pp=gcp;
    end

    t=-NVA.conv_t:NVA.conv_dt:NVA.conv_t; %ps
    posInds=t>0;
    
    tkMom=temporalKthMoment(rs, rd, k, optProp);
    RC=continuousReflectance(rs, rd, optProp);
    
    if NVA.usePar
        FFTconv=NVA.FFTconv;
        dt=NVA.conv_dt;
        l=NaN(size(r, 1), 1);
        parfor i=1:size(r, 1)
            lC=continuousPartPathLen(rs, r(i, :), rd, Vol, optProp);
            
            PHIsi=zeros(size(t));
            Rid=zeros(size(t));
            PHIsi(posInds)=temporalFluence(rs, r(i, :), t(posInds), ...
                optProp);
            Rid(posInds)=temporalReflectance(r(i, :), rd, t(posInds), ...
                optProp);
            
            if FFTconv
                convPHIsiRid=ifft(fft(PHIsi).*fft(Rid));
                convPHIsiRid((sum(posInds)+1):end)=[];
                
                l(i)=-(lC*tkMom-(Vol/RC)*...
                    sum(t(posInds).^k.*...
                    (convPHIsiRid*dt),...
                    2))*dt/tkMom;
            else
                l(i)=-(lC*tkMom-(Vol/RC)*...
                    trapz(t, t.^k.*...
                    (conv(PHIsi, Rid, 'same')*dt),...
                    2))/tkMom;
            end
        end
    else
        PHIsi=zeros(size(r, 1), length(t));
        Rid=zeros(size(r, 1), length(t));
        lC=continuousPartPathLen(rs, r, rd, Vol, optProp);
        PHIsi(:, posInds)=temporalFluence(rs, r, t(posInds), optProp);
        Rid(:, posInds)=temporalReflectance(r, rd, t(posInds), optProp);
        convPR=NaN(size(r, 1), length(t));
        convPHIsiRid=NaN(size(r, 1), sum(posInds));
        
        if NVA.FFTconv
            for i=1:size(r, 1)
                tmp=ifft(fft(PHIsi(i, :)).*fft(Rid(i, :)));
                convPHIsiRid(i, :)=tmp(1:sum(posInds));
            end
            
            dkMomIdmua=lC*tkMom-(Vol/RC)*...
                sum(...
                t(posInds).^k.*(convPHIsiRid*NVA.conv_dt)*NVA.conv_dt, 2);
        else
            for i=1:size(r, 1)
                convPR(i, :)=conv( ...
                    PHIsi(i, :), Rid(i, :), 'same')*NVA.conv_dt;
            end
            
            dkMomIdmua=lC*tkMom-(Vol/RC)*...
                trapz(t, ...
                t.^k.*convPR, 2);
        end
        
        l=-dkMomIdmua/tkMom;
    end
end