function S = FDTR_GetSHankel(kvect,xo,ws,A)
%definitions 
% kvect: matrix of wavevectors for evalution (NFreq columns do it for
% multiple frequencies at once)
% xo: beam offset distance (if xo > 4*ws, program will fail), (m)
% ws: spot size, 1/e^2 radius, (m)
% A: total beam intensity (Watts)  

%kvect is Nk x NFreq matrix
k=kvect(:,1);
kbar = k*ws/sqrt(2);
xbar = sqrt(2)*xo/ws;
if nargin ==3
    A =1;
end

prefactor = A/pi*exp(-(xbar^2+pi^2*kbar.^2));
x2 = xbar^2;
p = FDTR_GetPolyCoefs; %<---function that calculates polynomial coefficients based (Feser et al, RSI, 2012)
Nmax = length(p(:,1))-1;
Nk = length(kbar);

summand = zeros(size(k));
sigma = zeros(size(k));

for n = 0:Nmax
    PP=polyval(p(n+1,:),kbar);
    summand = x2^n/factorial(n)^2*PP;
    sigma = sigma + summand;
end

S = (prefactor.*sigma)*ones(1,length(kvect(1,:)));

    