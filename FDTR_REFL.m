%Calculates the Reflectivity Signal and Ratio
%In order to speed up the code, it is parallelized...the convention is...
%tdelay is COLUMN vector of desired delay times
%Mvect (the fourier components) are ROW vectors
%Matrices have size, length(tdelay) x length(Mvect)

%PARAMETERS (scalars unless told otherwise):  
%f:  modulation frequency 
%lambda:  VECTOR of thermal conductivities (W/mK) (layer 1 is the topmost layer)
%C:  VECTOR of specific heats (J/m3-K)
%t:  VECTOR of layer thicknesses (last layer is alway treated semi-inf, but
%       you still have to enter something).  (m)
%eta:  VECTOR of anisotropic ratio (kx/ky), use ones(length(lambda)) for
%       isotropic
%r_pump:  pump 1/e2 radius (m)
%r_probe: probe 1/e2 radius (m)

function [deltaT,phase]=FDTR_REFL(f,lambda,C,t,lambda_r,r_pump,r_probe,nnodes,lambda_anisotropy,lambda_ratio,xoffset)

% dT1=zeros(1,length(mvect))';
% dT2=zeros(1,length(mvect))';

kmax=1/sqrt(r_pump^2+r_probe^2)*4*pi;

%use Legendre-Gauss Integration
%computes weights and node locations...
if nargin<12
    nnodes = 35;
end
[kvect,weights]=lgwt_V4(nnodes,0,kmax);
% I1 = FDTR_TEMP(kvect,f,lambda,C,t,eta,r_pump,r_probe);
% I1 = FDTR_TEMP(kvect,f,lambda,C,t,eta,r_pump,r_probe,A_pump,xoffset);
I1=FDTR_TEMP(kvect,f,lambda,C,t,lambda_r,r_pump,r_probe,lambda_anisotropy,lambda_ratio,xoffset);
deltaT = weights'*I1;

phase=atan(imag(deltaT)./real(deltaT)).*180./pi;