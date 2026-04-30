function [R] = R_sphere_withPreCom(mu, preCom, NVA)
% [R] = R_sphere_withPreCom(mu, preCom, NVA)
% Angelo Sassaroli Ph.D. Fall 2023
% Giles Blaney Ph.D. Fall 2023

%calculation of the FD solution of DE in spherical geometry. AS October
%2023. For the geometry see A.Sassaroli et al. Appl. Opt. 40, 185-197 (2001).
%Here the Fourier Transform of Eq.(2) in the paper is calculated directly.
%Note that the unstable behavior that occurred in TD at times less than
%ballistic, here is reflected in the inaccuracy of AC and phase if the FT
%of the TPSF is carried out in the time range (0, Inf). So even in FD the
%calculation are correct only if the FT of the TPSF is carried out in the time range
%(tb, Inf), where tb is the ballistic time.
    
    arguments
        mu (1,2) double; %1/mm
        preCom struct;
    
        NVA.fmod (1,1) double = 140.625e6; %Hz
        NVA.ni (1,1) double = 1.4;
        NVA.no (1,1) double = 1;
    end
    c=0.299792458; %mm/ps
    
    mua=mu(1);
    musp=mu(2);
    
    rhos_arc=preCom.rhos_arc;
    rad=preCom.rad;
    
    nbf=preCom.nbf;
    nz=preCom.nz;
    
    z=preCom.bf_z;
    plg=preCom.plg;
    bsd_bn_05=preCom.bsd_bn_05;
    
    %% Calc Params
    theta=rhos_arc/rad; %rad
    D=1/(3*musp);
    v=c/NVA.ni;
    gsq=D*v;
    om=2*pi*NVA.fmod*10^-12; %rad/ps
    tb(1, 1, :)=2*rad*sin(theta./2)/v; %Ballistic time (ps)
    A=n2A(NVA.ni, NVA.no);
    
    dAA=A;
    z0=1/musp; 
    rs=rad-z0; %Radius of isotropic source
    re=rad+2*dAA*D; %Extrapolated radius
    
    bn_05=z./re; %beta_{n+1/2}
    kn=mua*v+gsq.*bn_05.^2-1i*om;
    
    %% Calc Complex Reflectance
    n_arr=repmat(1/2:1:nbf, nz, 1); %Bessel function index
    
    bs_bn_rs=besselj(n_arr, bn_05*rs);
    bs_bn_r=besselj(n_arr, bn_05*rad);
    bsd_bn_r=dbesselj(n_arr, bn_05*rad);
    
    S1=1./kn.*exp(-kn.*tb) ...
        .*bs_bn_rs./(bsd_bn_05.^2) ...
        .*(1/(2*rad)*bs_bn_r-bn_05.*bsd_bn_r); %FT (tb, Inf)
    S1=S1.*(1:2:(2*nbf)).*plg;
    S1=S1*gsq/(2*pi*re^2*sqrt(rad*rs));
    
    R=squeeze(sum(sum(S1)));
end

%% Supporting Functions
function dJndx = dbesselj(n,x)
%Derivative of Bessel function.
%It is based on the property dJn/dx=n/x*Jn(x)-J_{n+1}(x)
%which is equivalent to the property: dJn/dx=1/2*[J_{n-1}(x)-J_{n+1}(x)]
%
% DBESSELJ      A function that will generically calculate the
%               the derivative of a Bessel function of the first
%               kind of order n for all values of x.
%
% Example usage: dJndx = dbesselj(n,x);
% 
% INPUT ARGUMENTS
% ================
% n             Order of the Bessel function of the first kind
% x             Input variable to Bessel function
% 
% OUTPUT ARGUMENTS
% ================
% dJndx         Derivative of nth order Bessel function of the first 
%               kind at all values of x
    dJndx = n.*besselj(n, x)./x - besselj(n+1, x);
end