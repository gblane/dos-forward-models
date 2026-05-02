function [R] = R_FD_forward(rho, mua, musp, omega, v)
% [R] = R_FD_forward(rho, mua, musp, omega)
% Giles Blaney Spring 2019
% Inputs:
%   rho     - Source detector distance. (mm)
%   musp    - Reduced scattering. (1/mm)
%   mua     - Absorption. (1/mm)
%   omega   - Angular modulation frequency. (rad/sec)
%   v       - Speed of light in medium. (mm/sec)
%                
% Outputs:
%   R       - Complex reflectance. (1/mm^2)
    
    if nargin<=3
        omega=140.625e6*2*pi; %rad/sec
    end
    if nargin<=4
        c=2.99792458e11; %mm/sec
        v=c/1.4;
    end
    
    R=NaN(length(rho), 1);
    for i=1:length(rho)
        x0=1/musp; %mm
        rs=[x0, 0, 0];
        rd=[0, rho(i), 0];

        c=2.99792458e11; %mm/sec
        nin=c/v;
        nout=1;

        A=n2A(nin, nout);
%         D=1/(3*musp); %mm
%         xb=-2*A*D; %mm
        xb=-2*A/(3*(musp+mua)); %mm

%         mueff=sqrt(mua/D-1i*omega/(v*D)); %1/mm
        mueff=sqrt(3*(mua-1i*omega/v)*(musp+mua));
        
        rsp=[-x0+2*xb, 0, 0]; %mm

        r1=vecnorm(rd-rs, 2, 2);
        r2=vecnorm(rd-rsp, 2, 2);

        R(i)=(x0.*(1./r1+mueff).*exp(-mueff.*r1)./(r1.^2)+...
            (x0-2*xb).*(1./r2+mueff).*exp(-mueff.*r2)./(r2.^2))/...
            (4*pi);
    end

% % Equation 12.34 in Bigio and Fantini
% 
%     z0=1/musp;
%     
%     X=sqrt(3*(musp+mua))*sqrt(mua-1i*omega/c);
%     
%     R=(z0/(2*pi))*X*...
%         (1+1./(rho*X)).*...
%         exp(-rho*X)./(rho.^2);

end

