function [PHI] = complexFluence2L(rs, r, thk, en, optProp, opts, prtNm)
% complexFluence2L Calculate complex fluence in a two-layer medium.
%
% [PHI] = complexFluence2L(rs, r, thk, en, optProp, opts, prtNm)
%
% Written by Giles Blaney (Spring 2020; Ph.D. awarded May 2022)
%
% Inputs:
%   rs      - Source coordinates [mm]
%   r       - Coordinates to find fluence at [mm]
%   thk     - Layer thickness [mm]
%   en      - Bessel function roots [unitless]
%   optProp - Struct of optical properties [struct]
%   opts    - Options structure [struct]
%   prtNm   - Print name [string]
%
% Outputs:
%   PHI - Complex fluence [1/mm^2]
    
    tic;
    if nargin<=6
        prtNm='PHI';
    end
    
%     fprintf('Starting %s\n', prtNm);
    
    if nargin<=3
        optProp.nin=[1.4, 1.4];
        optProp.nout=1;
        optProp.musp=[1.20, 0.25]; %1/mm
        optProp.mua=[0.008, 0.020]; %1/mm
        
        opts.fmod=1.40625e8; %Hz
        opts.h_end=2000;
        opts.B=150; %mm
        en=zeroOrdBesselRoots(opts.h_end);
        
        warning(['Default optical properties used, this may be inconsistent'...
            ' with the musp used for source depth']);
    end
    
    if size(rs, 1)>1
        error('Can not use multiple sources');
    end
    
    % Make cylindrical geometry
    r_cyl=r-[rs(:, 1:2), zeros(size(rs(:, 3)))];
    rs_cyl=rs-[rs(:, 1:2), zeros(size(rs(:, 3)))];
    rho=sqrt(r_cyl(:, 1).^2+r_cyl(:, 2).^2);
    z=r_cyl(:, 3);
    lay1Inds=z<=thk;
    
    if size(en, 2)~=1
        en=en.';
    end
    sn=en(1:opts.h_end)/opts.B;

    c=2.99792458e11; %mm/sec
    nu=c./optProp.nin;
    omega=2*pi*opts.fmod; %rad/sec

    n=optProp.nin;
    A=n2A(n(1), optProp.nout);
    D=1./(3*optProp.musp); %mm
    zb=-2*A*D(1); %mm
    z0=rs_cyl(:, 3); %mm
    
    % [layer, Bessel zero]
    alpha=sqrt(optProp.mua./D+sn.^2+1i*omega./(D.*nu)).';
    
    G12_comden=(D(1)*alpha(1, :)*n(1)^2.*cosh(alpha(1, :)*(thk-zb))+...
        D(2)*alpha(2, :)*n(2)^2.*sinh(alpha(1, :)*(thk-zb)));
    G1_nonzdep=(sinh(alpha(1, :)*(z0-zb)).*...
        (n(1)^2*D(1)*alpha(1, :)-n(2)^2*D(2)*alpha(2, :)))./...
        ((D(1)*alpha(1, :).*exp(alpha(1, :)*(thk-zb))).*...
        G12_comden);
    G2_nonzdep=(n(2)^2*sinh(alpha(1, :)*(z0-zb)))./...
        G12_comden;
    J1sqr=(besselj(1, en(1:opts.h_end)).^2).';
    clear G12_comden;
    
    if size(r_cyl, 1)<=opts.h_end % CPU Limited
        PHI=nansum(...
            (lay1Inds.*(... G1 start
            ((exp(-alpha(1, :).*abs(z-z0))-exp(-alpha(1, :).*(z+z0-2*zb)))./...
            (2*D(1)*alpha(1, :)))+... G1 Term 1
            (sinh(alpha(1, :).*(z-zb)).*... G1 Term 2
            G1_nonzdep))+... G1 end
            ~lay1Inds.*(... G2 start
            exp(alpha(2, :).*(thk-z)).*...
            G2_nonzdep)).*... G2 end
            besselj(0, sn.'.*rho)./J1sqr,...
            2);
        PHI(isnan(PHI))=0;
    else % Memory Limited
        fprintf('Starting Bessel zeros loop for %s\n', prtNm);
        PHI=zeros(size(r_cyl, 1), 1);
        for i=1:opts.h_end
            PHItemp=...
                (lay1Inds.*(... G1 start
                ((exp(-alpha(1, i)*abs(z-z0))-exp(-alpha(1, i)*(z+z0-2*zb)))/...
                (2*D(1)*alpha(1, i)))+... G1 Term 1
                (sinh(alpha(1, i)*(z-zb))*... G1 Term 2
                G1_nonzdep(i)))+... G1 end
                ~lay1Inds.*(... G2 start
                exp(alpha(2, i)*(thk-z))*...
                G2_nonzdep(i))).*... G2 end
                besselj(0, sn(i)*rho)/J1sqr(i);
            PHItemp(isnan(PHItemp))=0;
            PHI=PHI+PHItemp;
            if mod(i, 100)==0
                fprintf('\t%s\t%d/%d\t%.3f sec\n',...
                    prtNm, i, opts.h_end, toc);
            end
        end
    end
    PHI=conj(PHI/(pi*opts.B^2));
    
%     fprintf('%s done %.3f sec\n', prtNm, toc);
end
