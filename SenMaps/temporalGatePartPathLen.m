function [l] = temporalGatePartPathLen(rs, r, rd, Vol, tg, optProp, NVA)
% Giles Blaney Ph.D. Spring 2023
% [l] = temporalGatePartPathLen(rs, rd, Vol, tg, optProp, NVA)
% Inputs:
%   rs      - Source corrdinates. (mm)
%   r       - Center corrdinate of volume. (mm)
%   rd      - Detector corrdinates. (mm)
%   Vol     - Volume. (mm^3)
%   tg      - Gate start and end time. (ps)
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
%   l       - Partial pathlength of gated t. (mm)
    
    arguments
        rs (1,3) double; %mm
        r (:,3) double; %mm
        rd (1,3) double; %mm
        Vol (1,1) double; %mm^3
        tg (2,:) double; %ps

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
    
    Rsd_t=temporalReflectance(rs, rd, t, optProp);
    
    tgInd=NaN(size(tg));
    Rsd_g=NaN(1, size(tg, 2));
    for j=1:size(tg, 2)
        [~, tgInd(1, j)]=min(abs(t-tg(1, j)));
        [~, tgInd(2, j)]=min(abs(t-tg(2, j)));
        i1=tgInd(1, j);
        i2=tgInd(2, j);
        
        Rsd_g(1, j)=trapz(t(i1:i2), Rsd_t(i1:i2), 2);
    end
    
    if NVA.usePar
        FFTconv=NVA.FFTconv;
        dt=NVA.conv_dt;
        l=NaN(size(r, 1), size(tg, 2));
        parfor i=1:size(r, 1)     
            
            PHIsi=zeros(size(t));
            Rid=zeros(size(t));
            PHIsi(posInds)=temporalFluence(rs, r(i, :), t(posInds), ...
                optProp);
            Rid(posInds)=temporalReflectance(r(i, :), rd, t(posInds), ...
                optProp);
            
            num=NaN(1, size(tg, 2));
            for j=1:size(tg, 2)
                i1=tgInd(1, j);
                i2=tgInd(2, j);
                
                if FFTconv
                    convPHIsiRid=ifft(fft(PHIsi).*fft(Rid));
                    convPHIsiRid((sum(posInds)+1):end)=[];
                    
                    [~, j1]=min(abs(t(i1)-t(posInds)));
                    [~, j2]=min(abs(t(i2)-t(posInds)));
                    
                    num(1, j)=sum(convPHIsiRid(j1:j2), 2)*dt;
                else
                    convTmp=conv(PHIsi, Rid, 'same')*dt;
                    num(1, j)=trapz(t(i1:i2), convTmp(i1:i2), 2);
                end
            end
            l(i, :)=(num./Rsd_g)*Vol;
        end
    else
        PHIsi=zeros(size(r, 1), length(t));
        Rid=zeros(size(r, 1), length(t));
        PHIsi(:, posInds)=temporalFluence(rs, r, t(posInds), optProp);
        Rid(:, posInds)=temporalReflectance(r, rd, t(posInds), optProp);
        convPR=NaN(size(r, 1), length(t));
        convPHIsiRid=NaN(size(r, 1), sum(posInds));
        
        if NVA.FFTconv
            for i=1:size(r, 1)
                tmp=ifft(fft(PHIsi(i, :)).*fft(Rid(i, :)));
                convPHIsiRid(i, :)=tmp(1:sum(posInds));
            end
            
            num=NaN(size(r, 1), size(tg, 2));
            for j=1:size(tg, 2)
                i1=tgInd(1, j);
                i2=tgInd(2, j);

                [~, j1]=min(abs(t(i1)-t(posInds)));
                [~, j2]=min(abs(t(i2)-t(posInds)));
    
                num(:, j)=sum(convPHIsiRid(:, j1:j2), 2)*NVA.conv_dt;
            end
        else
            for i=1:size(r, 1)
                convPR(i, :)=conv(PHIsi(i, :), Rid(i, :), 'same')*NVA.conv_dt;
            end
            
            num=NaN(size(r, 1), size(tg, 2));
            for j=1:size(tg, 2)
                i1=tgInd(1, j);
                i2=tgInd(2, j);
    
                num(:, j)=trapz(t(i1:i2), convPR(:, i1:i2), 2);
            end
        end
        l=(num./Rsd_g)*Vol;
    end
end