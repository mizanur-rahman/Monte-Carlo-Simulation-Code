s=[9.802 9.774;9.6305e-02 9.6305e-02];
se=[.037 .019;1.3805e-02 1.3805e-02 ];
for i=1:length(s)
    w(:,i)=1./(se(:,i).*se(:,i));
  
end
n=w.*s;
nu=sum(n');
den=sum(w');
mean=nu./den;
error=1./sqrt(den);