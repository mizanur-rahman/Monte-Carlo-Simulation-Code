function p = FDTR_GetPolyCoefs

Nmax = 40; %number of polynomials
Kmax = 2*Nmax+2;

p=zeros(Nmax,Kmax);
p(1,end)=pi;
L = [pi^2 0 -1 0];
M = [-1 0 1/(4*pi^2)];
N = [1/(4*pi^2) 0];
G = [-1 0];

for n=1:40 %normally 40
    index=n+1;
    pp = p(index-1,:);
    der = polyder(pp);
    secder = polyder(der);
    term1 = conv(L,pp);
    term2 = conv(M,der);
    term3 = conv(N,secder);
    NUM = polyadd(polyadd(term1,term2),term3);
   
    DEN = G;
     pnew = deconv(NUM,DEN);
    p(index,1:end)=pnew(end-Kmax+1:end);
end

end