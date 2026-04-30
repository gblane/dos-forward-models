function [R] = temporalReflectance(rs, rd, t, optProp, NVA)
% Giles Blaney Ph.D. Spring 2023
% [R] = temporalReflectance(rs, rd, tns, optProp)
% Inputs:
%   rs      - Isotropic source corrdinates. (mm)
%   rd      - Detector corrdinates. (mm)
%   t       - Time. (ps)
%   optProp - (OPTIONAL) Struct of optical properties with the following
%             fields:
%                nin  - (default=1.4) Index of refraction inside. (-)
%                nout - (default=1) Index of refraction outside. (-)
%                musp - (default=1.2 1/mm) Reduced scattering. (1/mm)
%                mua  - (default=0.01 1/mm) Absorption. (1/mm)
%   Name Value Arguments:
%           - 'simTyp' (default='DT'): String to switch between 'DT' and
%               'MC' simulation type.
%           - 'np' (default=1e9): Number of photons for MC.
%           - 'L' (default=[200, 200, 100]): Volume size for MC. (mm)
%           - 'grdSp' (default=100): Voxel side length for MC. (mm)
%           - 'g' (default=0.9): Anisotropy parameter.
%           - 'detRad' (default=1.5): Detector radius for MC. (mm)
% Outputs:
%   R       - Temporal reflectance. (1/(ps mm^2))
    
    arguments
        rs (:,3) double; %mm
        rd (:,3) double; %mm
        t (1,:) double; %ps

        optProp struct = [];

        NVA.simTyp (1,:) string = 'DT';

        % Only used in MC
        NVA.np (1,1) = 1e9;
        NVA.L (1,3) = [200, 200, 100]; %mm
        NVA.grdSp (1,1) = 50; %mm
        NVA.g (1,1) = 0.9;
        NVA.detRad (1,1) = 1.5; %mm
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

    nin=optProp.nin;
    nout=optProp.nout;
    musp=optProp.musp; %1/mm
    mua=optProp.mua; %1/mm
    
    posInds=t>0;
    tPos=t(posInds);
    
    c=0.299792458; %mm/ps
    v=c/nin;

    switch NVA.simTyp
        case 'DT'
            %% Run DT
            x0=rs(:, 1); %mm
            y0=rs(:, 2); %mm
            z0=rs(:, 3); %mm
            
            A=n2A(nin, nout);
            D=1/(3*musp); %mm
            zb=-2*A*D; %mm
        
            rsp=[x0, y0, -z0+2*zb]; %mm
        
            r1=vecnorm(rd-rs, 2, 2);
            r2=vecnorm(rd-rsp, 2, 2);
            
            R=(1/2)*...
                (exp(-mua*v*tPos)./((4*pi*D*v).^(3/2).*tPos.^(5/2))).*...
                ((z0.*exp(-r1.^2./(4*D*v*tPos)))+...
                ((z0-2*zb).*exp(-r2.^2./(4*D*v*tPos)))); %1/(ps mm^2)
        case 'MC'
            if size(rs, 1)>1
                error('Use only one source for MC');
            end

            %% Make MCX cfg0
            cfg0.seed=randi(255);
            
            % Req
            cfg0.nphoton=NVA.np;
            cfg0.maxdetphoton=1e8;
            cfg0.vol=uint8(ones(NVA.L/NVA.grdSp));
            cfg0.prop=[...
                0 0 1 nout;... Outside
                0 musp/(1-NVA.g) NVA.g nin... Inside
                ]; %mua mus g n
            
            tPos_sec=tPos*1e-12; %sec
            
            cfg0.tstart=median(diff(tPos_sec)); %sec
            cfg0.tstep=median(diff(tPos_sec));
            cfg0.tend=max(tPos_sec);
            
            xmid_src=size(cfg0.vol, 1)/2+1;
            ymid_src=size(cfg0.vol, 2)/2+1;
            
            cfg0.srcpos=[rs(1)/NVA.grdSp+xmid_src,...
                rs(2)/NVA.grdSp+ymid_src, 1];
            cfg0.srcdir=[0, 0, 1]; %vec
            cfg0.voidtime=0;
            
            % Detectors
            nDet=size(rd, 1);
            cfg0.detpos=[...
                rd(:, 1)/NVA.grdSp+xmid_src,...
                rd(:, 2)/NVA.grdSp+ymid_src,...
                ones(nDet, 1),...
                ones(nDet, 1)*NVA.detRad/NVA.grdSp];

            %Opt MC Sim
            cfg0.isreflect=1; %Consider n mismatch
            cfg0.unitinmm=NVA.grdSp; %mm
            
            %Opt GPU
            cfg0.autopilot=1;
            cfg0.gpuid=1;
            cfg0.isgpuinfo=1;
            
            %Opt SD
            cfg0.srctype='pencil';
            
            %Opt Output
            cfg0.issaveref=0;
            cfg0.outputtype='flux';
            cfg0.savedetflag='dpx';
            
            %Opt Debug
            cfg0.debuglevel='p';

            % Run MC
            [MCout, MCdetOut]=mcxlab(cfg0);
            if size(MCdetOut.ppath, 1)==cfg0.maxdetphoton
                warning('Max det photon reached, results may be incorrect');
            end

            % Post-process MC
            tMC=(cfg0.tstart:cfg0.tstep:cfg0.tend)*1e12; %ps
            
            E=MCout.stat.energytot;
            A=pi*NVA.detRad.^2;
            
            TPSF=NaN(nDet, length(tMC));
            R=NaN(nDet, length(tMC));
            for i=1:nDet
                detInds=MCdetOut.detid==i;
                
                l=MCdetOut.ppath(detInds)*NVA.grdSp;
                
                L=sum(l, 2);
                W=(E/cfg0.nphoton)*exp(-mua*L);

                Lt=L/v;

                for j=1:length(tMC)
                    inds=and(Lt>tMC(j), Lt<=(tMC(j)+cfg0.tstep*1e12));
                    TPSF(i, j)=sum(W(inds));
                end
                R(i, :)=TPSF/(E*A*cfg0.tstep*1e12); %1/(ps mm^2)
            end

        otherwise
            error('Unkown simTyp');
    end
    R=[zeros(size(R, 1), sum(~posInds)), R];
end