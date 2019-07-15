%Computes frequency domain average temperature response to 
%Periodic gaussian pump beam, probed by another gaussian beam
%
%This program is vectorized/optimized to handle all frequencies 
%simultaneously (f is a ROW vector)
%
%Definitions
%f: excitation frequency (Hz), ROW vector
%lambda: vector of thermal conductivities, 
%   lambda(1)=top surface,(W/m-K)
%C: vector of volumetric specific heat (J/m3-K)
%t: thicknesses of each layer (layer N will NOT be used, semiinfinite)
%r_pump: Pump spot size (m)
%r_probe: Probe spot size (m)
%A_pump: Pump power (W), used to ESTIMATE amplitude (not used for fitting)

function Integrand=FDTR_TEMP(kvectin,freq,lambda,C,t,lambda_r,r_pump,r_probe,lambda_anisotropy,lambda_ratio,xoffset)
A_pump=1;
Nfreq=length(freq);
kvect=kvectin(:)*ones(1,Nfreq);
Nlayers=length(lambda); %# of layers
Nint=length(kvectin); %# of different frequencies to calculate for

%k is a COLUMN vector (actually a matrix that changes down the rows)
%f is a ROW vector

%kmax=1/sqrt(r_pump^2+r_probe^2)*1.5; %cutoff wavevector for integration
ii=sqrt(-1);
alpha=lambda./C;
omega=2*pi*freq;
q2=ones(Nint,1)*(ii*omega./alpha(Nlayers));
kvect2=kvect.^2;


un=sqrt(4*pi^2*(lambda_r(Nlayers)./lambda(Nlayers))*kvect2+q2);
gamman=lambda(Nlayers)*un;
Bplus=zeros(Nint,Nfreq);
Bminus=ones(Nint,Nfreq);
kterm2=4*pi^2*kvect2;
if Nlayers~=1
    for n=Nlayers:-1:2
        q2=ones(Nint,1)*(ii*omega./alpha(n-1));
        unminus=sqrt((lambda_r(n-1)./lambda(n-1))*kterm2+q2);
        gammanminus=lambda(n-1)*unminus;
        AA=gammanminus+gamman;
        BB=gammanminus-gamman;
        temp1=AA.*Bplus+BB.*Bminus;
        temp2=BB.*Bplus+AA.*Bminus;
        expterm=exp(unminus*t(n-1));
        Bplus=(0.5./(gammanminus.*expterm)).*temp1;
        Bminus=0.5./(gammanminus).*expterm.*temp2;
        % These next 3 lines fix a numerical stability issue if one of the
        % layers is very thick or resistive;
        penetration_logic=logical(t(n-1)*abs(unminus)>100);  %if pentration is smaller than layer...set to semi-inf
        Bplus(penetration_logic)=0;
        Bminus(penetration_logic)=1;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        un=unminus;
        gamman=gammanminus;
    end
end

G=(Bplus+Bminus)./(Bminus-Bplus)./gamman; %The layer G(k)
P=A_pump*exp(-pi^2*r_pump^2/2*kvect.^2);
S = FDTR_GetSHankel(kvect,xoffset,r_probe,1);

Integrand=G.*P.*S.*kvect; %The rest of the integrand

        
