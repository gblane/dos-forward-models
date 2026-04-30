function [preCom] = get_Rshpere_preCom(rhos_arc, rad, NVA)
% [preCom] = get_Rshpere_preCom(rho_arc, rad)
% Supporting function for R_sphere_withPreCom()
% Giles Blaney Ph.D. Fall 2023
    
    arguments
        rhos_arc (1,:) double; %Source-Detector Arc Lengths (mm)
        rad (1,1) double; %Sphere Radius (mm)

        NVA.nbf (1,1) double = 200; %Number of Bessel functions
        NVA.nz (1,1) double = 200; %Number of zeros of each Bessel function
    end
    
    theta(1, 1, :)=rhos_arc/rad; %rad
    
    %% Bessel Zeros
    n_arr=(1/2:1:NVA.nbf)'; %Bessel function index
    bf_z=besselzero(n_arr, NVA.nz, 1)';
    
    %% rho_arc Dependent Values
    inds=repmat(0:(NVA.nbf-1), 1, 1, size(theta, 3));
    cosTh=repmat(cos(theta), 1, NVA.nbf, 1);
    plg=legendreP(inds, cosTh);
    
    n_arr=repmat(1/2:1:NVA.nbf, NVA.nz, 1);
    bsd_bn_05=dbesselj(n_arr, bf_z);
    
    %% Package
    preCom.nbf=NVA.nbf;
    preCom.nz=NVA.nz;
    
    preCom.bf_z=bf_z; %Bessel zeros
    preCom.plg=plg;
    preCom.bsd_bn_05=bsd_bn_05;
    
    preCom.rhos_arc=rhos_arc;
    preCom.rad=rad;
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