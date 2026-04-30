function [S, params, Svox] = makeS(typStr, rs, rd, optProp, NVA)
% [S, params, Svox] = makeS(typStr, rs, rd, optProp, NVA)
% 
% Giles Blaney Ph.D. Summer 2023
% 
% Inputs:
%   typStr   - String setting the type of data considered in the format:
%               AA_BB_C where, [e.g., CW_SD_I or FD_DS_P]
%               AA is the temporal domain either: 
%                   CW  - Continuous Wave
%                   FD  - Frequency-Domain 
%                   TD  - Time-Domain 
%               BB is the arrangement type either:
%                   SD  - Single-Distance 
%                   SS  - Single-Slope 
%                   DS  - Dual-Slope 
%               CC is the data-type either:
%                   I   - Intensity for CW or FD (Amplitude for FD)
%                   GI  - Gated Intensity for TD
%                   DGI - Difference in Gated Intensity for TD 
%                   P   - Phase for FD
%                   T   - Mean time-of-flight for TD
%                   V   - Variance for TD
%   rs       - Pencil beam source corrdinates [x, y, z]. (mm)
%   rd       - Detector corrdinates in [x, y, z]. (mm)
%   optProp  - (OPTIONAL) Struct of optical properties with the following
%                fields:
%                    nin  - (default=1.333) Index of refraction inside. (-)
%                    nout - (default=1) Index of refraction outside. (-)
%                    musp - (default=1.1 1/mm) Reduced scattering. (1/mm)
%                    g    - (default=0.9) Anisotropy for MC. (-)
%                    mua  - (default=0.011 1/mm) Absorption. (1/mm)
% 
% Name Value Arguments:
%   'pert'   - (default=[1, 1, 1] mm) Size of perturbation in [x, y, z], 
%               must be a multiple of voxel size (i.e., dr). (mm)
%   'dr'     - (default=1 mm) Voxel side length. (mm)
%   'xl'     - (default: 10 mm past the extent of the optodes) Limits in the
%               x direction in [min, max]. (mm)
%   'yl'     - (default: 10 mm past the extent of the optodes) Limits in the
%               y direction in [min, max]. (mm)
%   'zl'     - (default=[0, 25] mm) Limits in the y direction in [min, max].
%               (mm)
%   'fmod'   - (default=100e6 Hz) Modulation frequency for FD. (Hz)
%   'ndt'    - (default=10e3) Number of time bins for TD or MC.
%   'tend'   - (default=10e3 ps) Max time bin edge for TD or MC. (ps) 
%   'tg'     - (default=[1000, 2000] ps) Gate start and end time for TD GI.
%               (ps)
%   'tgE     - (default=[0, 1000] ps) Early gate start and end time for 
%               TD DGI. (ps)
%   'simTyp' - (default='DT'): String to switch between 'DT' 
%               (Diffusion Theory) and 'MC' (Monte Carlo w/ MCX) 
%               simulation type.
%   'detNA'  - (default=0.5) Detector NA for MC.
%   'np'     - (default=1e8): Number of photons for MC.
%   'usePar' - (default=true): Use parfoor loops.
%   'FFTconv'- (default=true): Use ifft(fft*fft) as conv.
% 
% Outputs:
%   S        - Sensitivity (unit-less)
%   params   - Struct with the following fields
%               x        - X-axis vector. (mm)
%               y        - Y-axis vector. (mm)
%               z        - Z-axis vector. (mm)
%               xl_valid - Valid x-axis limits with perturbation size. (mm)
%               yl_valid - Valid y-axis limits with perturbation size. (mm)
%               zl_valid - Valid z-axis limits with perturbation size. (mm)
%   Svox     - Voxelized sensitivity not considering the perturbation size.
%               (unit-less)

    arguments
        typStr (1,:) string;
        rs (:,3) double; %mm [x, y, z]
        rd (:,3) double; %mm [x, y, z]
        
        optProp struct = [];
        
        NVA.pert (1,3) double = [1, 1, 1]; %mm
        NVA.dr (1,1) double = 1; %mm
        NVA.xl (1,2) double = [NaN, NaN]; %mm
        NVA.yl (1,2) double = [NaN, NaN]; %mm
        NVA.zl (1,2) double = [0, 25]; %mm
        
        NVA.fmod (1,1) double = 100e6; %Hz
        
        NVA.ndt (1,1) double = 10e3;
        NVA.tend (1,1) double = 10e3; %ps

        NVA.tg (2,1) double = [1000; 2000]; %ps
        NVA.tgE (2,1) double = [500; 1000]; %ps

        NVA.simTyp (1,:) string = 'DT';
        NVA.detNA (1,1) double = 0.5; 

        NVA.np (1,1) double = 1e8;

        NVA.usePar (1,1) logical = true;
        NVA.FFTconv (1,1) logical = true;
    end
    
    typ=split(upper(typStr), '_');
    
    %% Check Inputs
    % Default optical properties if needed
    if isempty(optProp)
        clear optProp;

        optProp.nin=1.333;
        optProp.nout=1;
        optProp.musp=1.1; %1/mm
        optProp.mua=0.011; %1/mm
        
        warning('Default optical properties used');
    end
    
    % Check rs and rd based on arrangement type
    % Then make iso source and rep for num meas pairs
    if strcmp(NVA.simTyp, 'MC')
        z0=0;
    else
        z0=1/optProp.musp;
    end
    switch typ{2}
        case 'SD'
            if size(rs, 1)~=1 || size(rd, 1)~=1
                error('Incorrect number of optodes for %s', typ{2});
            end
            rSrcs=rs+[0, 0, z0];
            rDets=rd;
        case 'SS'
            if (size(rs, 1)==1 && size(rd, 1)==2)
                rSrcs=[rs; rs]+[0, 0, z0];
                rDets=rd;
            elseif (size(rs, 1)==2 && size(rd, 1)==1)
                rSrcs=rs+[0, 0, z0];
                rDets=[rd; rd];
            else
                error('Incorrect number of optodes for %s', typ{2});
            end
            
        case 'DS'
            if size(rs, 1)~=2 || size(rd, 1)~=2
                error('Incorrect number of optodes for %s', typ{2});
            end
            rSrcs=[rs([1, 1], :); rs([2, 2], :)]+[0, 0, z0];
            rDets=[rd; flipud(rd)];
        otherwise
            error('Unknown arrangement type %s', typ{2});
    end
    
    % Check pert is a multiple of dr
    if any(mod(NVA.pert, NVA.dr)~=0)
        error('pert size must be a multiple of dr')
    end

    % Default xl and yl if needed
    if any(isnan(NVA.xl))
        NVA.xl=[min([rs(:, 1); rd(:, 1)])-10, ...
            max([rs(:, 1); rd(:, 1)])+10];
    end
    if any(isnan(NVA.yl))
        NVA.yl=[min([rs(:, 2); rd(:, 2)])-10, ...
            max([rs(:, 1); rd(:, 2)])+10];
    end
    
    % Check if MC can be run
    if strcmp(NVA.simTyp, 'MC')
        if strcmp(typ{1}, 'FD')
            error('MC not supported for FD at this time');
        end

        if ~(strcmp(typ{3}, 'I') || strcmp(typ{3}, 'GI') || ...
                strcmp(typ{3}, 'DGI'))
            error('MC only supported for I and GI at this time');
        end
    end
    
    %% Check if run MC and make coordinates
    if strcmp(NVA.simTyp, 'MC')
        for i=1:size(rSrcs, 1)
            [adjoint(i), MCXparams(i), ~]=myMCXLAB_adjoint( ...
                rSrcs(i, :), rDets(i, :), optProp, ...
                'np', NVA.np, 'dr', NVA.dr, ...
                'xl', NVA.xl, 'yl', NVA.yl, 'zl', NVA.zl, ...
                'ndt', NVA.ndt, 'tend', NVA.tend, ...
                'detNA', NVA.detNA);
        end
        XX=MCXparams(1).XX;
        params.x=MCXparams(1).x;
        params.y=MCXparams(1).y;
        params.z=MCXparams(1).z;
    else        
        [YY, XX, ZZ]=meshgrid( ...
            NVA.yl(1):NVA.dr:NVA.yl(2), ...
            NVA.xl(1):NVA.dr:NVA.xl(2), ...
            NVA.zl(1):NVA.dr:NVA.zl(2));
        params.x=squeeze(XX(:, 1, 1));
        params.y=squeeze(YY(1, :, 1));
        params.z=squeeze(ZZ(1, 1, :));
        r=[XX(:), YY(:), ZZ(:)];
    end

    %% Switch over measurement type and calculate generalized path-lengths
    switch typ{1}
        case 'CW'
            switch typ{3}
                case 'I'
                    L=NaN(size(rSrcs, 1), 1);
                    Y=NaN(size(rSrcs, 1), 1);
                    ll=NaN([size(XX), size(rSrcs, 1)]);
                    for i=1:size(rSrcs, 1)
                        Y(i)=1;
                        switch NVA.simTyp
                            case 'DT'
                                L(i)=continuousTotPathLen(...
                                    rSrcs(i, :), rDets(i, :), optProp);
                                
                                l=continuousPartPathLen( ...
                                    rSrcs(i, :), r, rDets(i, :), ...
                                    NVA.dr^3, optProp);
                                l(isnan(l))=0;
                                ll(:, :, :, i)=reshape(l, size(XX));
                            case 'MC'
                                [ll(:, :, :, i), L(i)]=...
                                    continuousPathLen_MCadjoint(...
                                    adjoint(i), MCXparams(i).t, NVA.dr^3);
                            otherwise
                                error('Unknown simTyp %s', NVA.simTyp);
                        end
                    end
                otherwise
                    error('Unkown data-type %s for %s', typ{3}, typ{1});
            end
        
        case 'FD'
            switch typ{3}
                case 'I'
                    L=NaN(size(rSrcs, 1), 1);
                    Y=NaN(size(rSrcs, 1), 1);
                    ll=NaN([size(XX), size(rSrcs, 1)]);
                    for i=1:size(rSrcs, 1)
                        L(i)=real(complexTotPathLen(...
                            rSrcs(i, :), rDets(i, :), ...
                            2*pi*NVA.fmod, optProp));
                        Y(i)=1;
                        
                        l=real(complexPartPathLen( ...
                            rSrcs(i, :), r, rDets(i, :), ...
                            NVA.dr^3, 2*pi*NVA.fmod, optProp));
                        l(isnan(l))=0;
                        
                        ll(:, :, :, i)=reshape(l, size(XX));
                    end
                case 'P'                   
                    L=NaN(size(rSrcs, 1), 1);
                    Y=NaN(size(rSrcs, 1), 1);
                    ll=NaN([size(XX), size(rSrcs, 1)]);
                    for i=1:size(rSrcs, 1)
                        L(i)=imag(complexTotPathLen(...
                            rSrcs(i, :), rDets(i, :), ...
                            2*pi*NVA.fmod, optProp));
                        Y(i)=1;

                        l=imag(complexPartPathLen( ...
                            rSrcs(i, :), r, rDets(i, :), ...
                            NVA.dr^3, 2*pi*NVA.fmod, optProp));
                        l(isnan(l))=0;
                        
                        ll(:, :, :, i)=reshape(l, size(XX));
                    end
                otherwise
                    error('Unkown data-type %s for %s', typ{3}, typ{1});
            end
        
        case 'TD'
            switch typ{3}
                case 'GI'                    
                    L=NaN(size(rSrcs, 1), 1);
                    Y=NaN(size(rSrcs, 1), 1);
                    ll=NaN([size(XX), size(rSrcs, 1)]);
                    for i=1:size(rSrcs, 1)
                        Y(i)=1;
                        switch NVA.simTyp
                            case 'DT'
                                L(i)=temporalGateTotPathLen( ...
                                    rSrcs(i, :), rDets(i, :), ...
                                    NVA.tg, optProp, ...
                                    'conv_t', NVA.tend, ...
                                    'conv_dt', NVA.tend/NVA.ndt);
        
                                l=temporalGatePartPathLen( ...
                                    rSrcs(i, :), r, rDets(i, :), ...
                                    NVA.dr^3, NVA.tg, optProp, ...
                                    'conv_t', NVA.tend, ...
                                    'conv_dt', NVA.tend/NVA.ndt, ...
                                    'usePar', NVA.usePar, ...
                                    'FFTconv', NVA.FFTconv);
                                l(isnan(l))=0;
                                ll(:, :, :, i)=reshape(l, size(XX));
                            case 'MC'
                                [ll(:, :, :, i), L(i)]=...
                                    temporalGatePathLen_MCadjoint(...
                                    adjoint(i), MCXparams(i).t, NVA.dr^3,...
                                    NVA.tg);
                            otherwise
                                error('Unknown simTyp %s', NVA.simTyp);
                        end
                    end
                case 'DGI'
                    if size(rSrcs, 1)~=1
                        error('DGI only supported for SD at this time');
                    end
                    typ{2}=[typ{2}, '_DIFF'];
                    tg_all=[NVA.tgE, NVA.tg];

                    L=NaN(size(tg_all, 2), 1);
                    Y=NaN(size(tg_all, 2), 1);
                    ll=NaN([size(XX), size(tg_all, 2)]);
                    for i=1:size(tg_all, 2)
                        Y(i)=1;
                        switch NVA.simTyp
                            case 'DT'
                                L(i)=temporalGateTotPathLen(...
                                    rSrcs, rDets, ...
                                    tg_all(:, i), optProp, ...
                                    'conv_t', NVA.tend, ...
                                    'conv_dt', NVA.tend/NVA.ndt);
        
                                l=temporalGatePartPathLen( ...
                                    rSrcs, r, rDets, ...
                                    NVA.dr^3, tg_all(:, i), optProp, ...
                                    'conv_t', NVA.tend, ...
                                    'conv_dt', NVA.tend/NVA.ndt, ...
                                    'usePar', NVA.usePar, ...
                                    'FFTconv', NVA.FFTconv);
                                l(isnan(l))=0;
                                ll(:, :, :, i)=reshape(l, size(XX));
                            case 'MC'
                                [ll(:, :, :, i), L(i)]=...
                                    temporalGatePathLen_MCadjoint(...
                                    adjoint, MCXparams.t, NVA.dr^3,...
                                    tg_all(:, i));
                            otherwise
                                error('Unknown simTyp %s', NVA.simTyp);
                        end
                    end
                case 'T'                    
                    L=NaN(size(rSrcs, 1), 1);
                    Y=NaN(size(rSrcs, 1), 1);
                    ll=NaN([size(XX), size(rSrcs, 1)]);
                    for i=1:size(rSrcs, 1)
                        L(i)=temporalKthMomTotPathLen(...
                            rSrcs(i, :), rDets(i, :), 1, optProp);
                        Y(i)=temporalKthMoment(...
                            rSrcs(i, :), rDets(i, :), 1, optProp);
                        
                        l=temporalKthMomPartPathLen( ...
                            rSrcs(i, :), r, rDets(i, :), ...
                            NVA.dr^3, 1, optProp, ...
                            'conv_t', NVA.tend, ...
                            'conv_dt', NVA.tend/NVA.ndt, ...
                            'usePar', NVA.usePar, ...
                            'FFTconv', NVA.FFTconv);
                        l(isnan(l))=0;
                        
                        ll(:, :, :, i)=reshape(l, size(XX));
                    end
                case 'V'                    
                    L=NaN(size(rSrcs, 1), 1);
                    Y=NaN(size(rSrcs, 1), 1);
                    ll=NaN([size(XX), size(rSrcs, 1)]);
                    for i=1:size(rSrcs, 1)
                        L(i)=temporalVarTotPathLen(...
                            rSrcs(i, :), rDets(i, :), optProp);
                        Y(i)=temporalVar(...
                            rSrcs(i, :), rDets(i, :), optProp);
                        
                        l=temporalVarPartPathLen( ...
                            rSrcs(i, :), r, rDets(i, :), ...
                            NVA.dr^3, optProp, ...
                            'conv_t', NVA.tend, ...
                            'conv_dt', NVA.tend/NVA.ndt, ...
                            'usePar', NVA.usePar, ...
                            'FFTconv', NVA.FFTconv);
                        l(isnan(l))=0;
                        
                        ll(:, :, :, i)=reshape(l, size(XX));
                    end
                otherwise
                    error('Unkown data-type %s for %s', typ{3}, typ{1});
            end
            
        otherwise
            error('Unknown temporal domain %s', typ{1});
    end
    
    %% Switch over arrangement type and calculate S
    switch typ{2}
        case 'SD'
            Svox=ll/L;
        case {'SS', 'SD_DIFF'}
            Svox=(Y(2)*ll(:, :, :, 2)-Y(1)*ll(:, :, :, 1))/...
                (Y(2)*L(2)-Y(1)*L(1));
        case 'DS'
            Svox=((Y(2)*ll(:, :, :, 2)-Y(1)*ll(:, :, :, 1))+...
                (Y(4)*ll(:, :, :, 4)-Y(3)*ll(:, :, :, 3)))/...
                ((Y(2)*L(2)-Y(1)*L(1))+...
                (Y(4)*L(4)-Y(3)*L(3)));
        otherwise
            error('Unknown arrangement type %s', typ{2});
    end

    %% Apply pert size
    H=ones(NVA.pert/NVA.dr);
    S=convn(Svox, H, 'same');
    
    if all(NVA.pert==1)
        params.xl_valid=params.x([1, end]);
        params.yl_valid=params.y([1, end]);
        params.zl_valid=params.z([1, end]);
    else
        params.xl_valid=[1, -1]*NVA.pert(1)/2+params.x([1, end]);
        params.yl_valid=[1, -1]*NVA.pert(2)/2+params.y([1, end]);
        params.zl_valid=[1, -1]*NVA.pert(3)/2+params.z([1, end]);
    end

end